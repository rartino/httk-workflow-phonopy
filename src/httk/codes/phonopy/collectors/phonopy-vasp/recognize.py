"""Recognize hook for ``phonopy.calculation.vasp``: a phonopy displacement set and its runs."""

import hashlib
from pathlib import Path

from httk.workflow.calculations import content_digest
from httk.workflow.collecting import existing_file
from httk.workflow.hookapi import Claim, Unclaimed

from httk.codes.phonopy.collect import displacement_directories, read_displacement_count


def recognize(directory: Path) -> Claim | Unclaimed | None:
    """Claim a phonopy directory with one ``disp-NNN`` subdirectory per displacement.

    :param directory: The directory holding ``phonopy_disp.yaml`` and the unit-cell ``POSCAR``.
    :return: A claim that consumes the displacement directories, or ``Unclaimed``
        when the file cannot be read or the directories do not number the
        displacements ``1`` to ``N``. The claim's identity covers
        ``phonopy_disp.yaml`` and ``POSCAR`` and, per displacement directory, the
        VASP inputs that identify its run, so the same displacement set computed
        with other VASP settings is another calculation.
    """

    path = existing_file(directory / "phonopy_disp.yaml")
    if path is None:
        return None
    try:
        count = read_displacement_count(path)
    except ValueError as exc:
        return Unclaimed(str(exc))
    found = displacement_directories(directory)
    if not found or len(found) != count:
        return Unclaimed(f"{len(found)} displacement directories for {count} displacements")
    if [int(name.removeprefix("disp-")) for name in found] != list(range(1, count + 1)):
        return Unclaimed(f"displacement directories {', '.join(found)} are not numbered 1 to {count}")
    identity = hashlib.sha256(content_digest(directory, ["phonopy_disp.yaml", "POSCAR"]).encode() + b"\0")
    for name in found:
        # The input set the vasp.calculation.static collector identifies each run by.
        inputs = content_digest(directory / name, ["INCAR", "KPOINTS", "POTCAR", "POSCAR"])
        identity.update(name.encode() + b"\0" + inputs.encode() + b"\0")
    return Claim(identity.hexdigest(), consumes=found)
