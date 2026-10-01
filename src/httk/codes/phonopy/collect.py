"""Building blocks for collecting finite-displacement phonopy calculations.

phonopy's finite-displacement method writes the displaced supercells of a unit
cell (``phonopy -d``) and records the displacements in ``phonopy_disp.yaml``.
Each displaced supercell is then computed by a force calculator in its own
directory, conventionally ``disp-001``, ``disp-002``, and so on. These helpers
read the number of displacements and find those directories, for the
``phonopy.calculation.vasp`` collector and for collectors of one's own.
"""

import re
from pathlib import Path
from typing import IO, cast

from httk.core.datastream.compression import open_compressed

__all__ = ["displacement_directories", "read_displacement_count"]

_SECTION = re.compile(r"displacements:[ \t]*(#.*)?")
_ITEM = re.compile(r"-([ \t]|$)")
_ENTRY = re.compile(r"-[ \t]+(atom:[ \t]*([0-9]+)|#[ \t]*[0-9]+)[ \t]*(#.*)?")
_DIRECTORY = re.compile(r"disp-([0-9]+)")


def read_displacement_count(path: Path) -> int:
    """Count the displacements recorded in a phonopy ``phonopy_disp.yaml``.

    The file is scanned line by line rather than parsed as YAML, for the layouts
    phonopy itself writes: a top-level ``displacements:`` key with a block list
    of one entry per displaced supercell. An entry of phonopy's usual dataset
    (``phonopy -d``) starts with ``- atom: N``, the 1-based index of the
    displaced supercell atom, and continues with indented keys such as
    ``displacement:``. An entry of a random-displacement dataset
    (``phonopy -d --rd N``) starts with ``- # N`` and continues with an indented
    list of ``- displacement:`` items. The entries are the list items at the
    indentation of the first one (column 0 as phonopy writes them, or deeper as
    YAML rewriters such as ruamel.yaml indent them); deeper lines belong to an
    entry. The list ends at a line indented less than its entries or at the next
    top-level key. Other sections, comments and blank lines are skipped, and the
    file may be compressed (its name then carries the compression suffix, as in
    ``phonopy_disp.yaml.gz``).

    :param path: The ``phonopy_disp.yaml`` file, possibly compressed.
    :return: The number of displacements, which is zero for an empty list.
    :raises ValueError: If the file has no top-level ``displacements:`` block list,
        has more than one, or a list entry starts with neither ``- atom: N``
        (``N >= 1``) nor ``- # N``, or the file is not UTF-8 text.
    """

    count: int | None = None
    inside = False
    indent: int | None = None  # of the entries of the displacements list
    with path.open("rb") as raw, open_compressed(raw, compression="extension", name=path.name) as stream:
        for number, data in enumerate(cast(IO[bytes], stream), start=1):
            line = data.decode("utf-8-sig" if number == 1 else "utf-8").rstrip("\r\n").rstrip()
            content = line.lstrip(" \t")
            if not content or content.startswith("#"):
                continue
            depth = len(line) - len(content)
            item = _ITEM.match(content) is not None
            if inside:
                if indent is None and item:
                    indent = depth
                if indent is None and depth > 0:
                    raise ValueError(f"{path.name}:{number}: 'displacements:' is not a block list")
                if indent is not None and depth < indent:
                    inside = False
                elif item and depth == indent:
                    entry = _ENTRY.fullmatch(content)
                    if entry is None or (entry.group(2) is not None and int(entry.group(2)) < 1):
                        raise ValueError(
                            f"{path.name}:{number}: a displacements entry must start with '- atom: N' or '- # N': "
                            f"{content!r}"
                        )
                    count = cast(int, count) + 1
            if depth > 0 or item:
                continue
            # A top-level key: it ends any displacements list and may start one.
            inside = _SECTION.fullmatch(line) is not None
            if inside:
                if count is not None:
                    raise ValueError(f"{path.name}:{number}: a second top-level 'displacements:' list")
                count, indent = 0, None
    if count is None:
        raise ValueError(f"{path.name} has no top-level 'displacements:' list; is it a phonopy_disp.yaml?")
    return count


def displacement_directories(directory: Path) -> tuple[str, ...]:
    """Find the displacement subdirectories of a phonopy calculation directory.

    These are the subdirectories named ``disp-`` followed by decimal digits,
    phonopy's usual naming (``disp-001`` for the supercell of ``POSCAR-001``).
    They are ordered by their number, and by name among equal numbers.

    :param directory: The phonopy calculation directory.
    :return: The subdirectory names, in displacement order.
    """

    found: list[tuple[int, str]] = []
    for path in Path(directory).iterdir():
        match = _DIRECTORY.fullmatch(path.name)
        if match is not None and path.is_dir():
            found.append((int(match.group(1)), path.name))
    return tuple(name for _number, name in sorted(found))
