"""The building blocks of :mod:`httk.codes.phonopy.collect`."""

import bz2
import gzip
import lzma
import re
from collections.abc import Callable
from pathlib import Path
from typing import IO, Any

import pytest

from conftest import phonopy_disp_yaml
from httk.codes.phonopy.collect import displacement_directories, read_displacement_count


def _write(directory: Path, text: str, name: str = "phonopy_disp.yaml") -> Path:
    path = directory / name
    path.write_text(text, encoding="utf-8")
    return path


def test_the_displacements_of_a_phonopy_disp_yaml_are_counted(tmp_path: Path) -> None:
    assert read_displacement_count(_write(tmp_path, phonopy_disp_yaml())) == 2
    assert read_displacement_count(_write(tmp_path, phonopy_disp_yaml((1, 9, 2, 10, 3)))) == 5


def test_a_random_displacement_dataset_is_counted(tmp_path: Path) -> None:
    # phonopy -d --rd N: entries start with "- # N" and list every atom's displacement.
    assert read_displacement_count(_write(tmp_path, phonopy_disp_yaml(random=True))) == 2
    assert read_displacement_count(_write(tmp_path, phonopy_disp_yaml((1, 2, 3, 4, 5), random=True))) == 5


def test_an_indented_list_is_counted_at_the_indent_of_its_first_entry(tmp_path: Path) -> None:
    # As ruamel.yaml rewrites a file with a sequence indent of 2.
    text = phonopy_disp_yaml()
    head, _, tail = text.partition("displacements:\n")
    indented = head + "displacements:\n" + "".join(f"  {line}\n" for line in tail.splitlines())
    assert read_displacement_count(_write(tmp_path, indented + "later:\n- atom: 4\n")) == 2
    text = "displacements:\n    - atom: 1\n      displacement: [0.01, 0, 0]\n  - atom: 2\n"
    assert read_displacement_count(_write(tmp_path, text)) == 1


def test_any_blank_after_the_dash_is_accepted(tmp_path: Path) -> None:
    assert read_displacement_count(_write(tmp_path, "displacements:\n-\tatom: 1\n-   atom: 2\n-  #   3\n")) == 3


def test_indented_lists_of_an_entry_are_not_counted(tmp_path: Path) -> None:
    # phonopy_params.yaml carries each displacement's forces as an indented list.
    assert read_displacement_count(_write(tmp_path, phonopy_disp_yaml(forces=True))) == 2


@pytest.mark.parametrize(
    ("suffix", "opener"), [(".gz", gzip.open), (".bz2", bz2.open), (".xz", lzma.open)], ids=["gz", "bz2", "xz"]
)
def test_a_compressed_file_is_read(tmp_path: Path, suffix: str, opener: Callable[..., IO[Any]]) -> None:
    path = tmp_path / f"phonopy_disp.yaml{suffix}"
    with opener(path, "wb") as handle:
        handle.write(phonopy_disp_yaml().encode())
    assert read_displacement_count(path) == 2


def test_the_list_ends_at_the_next_top_level_key(tmp_path: Path) -> None:
    text = "displacements:\n# a comment\n\n- atom: 3 # three\n  displacement:\n    [ 0.01, 0, 0 ]\nlater:\n- atom: 4\n"
    assert read_displacement_count(_write(tmp_path, text)) == 1
    assert read_displacement_count(_write(tmp_path, "unit_cell:\n  lattice: []\ndisplacements:\nlater: 1\n")) == 0


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("phonopy:\n  version: 2\n", "no top-level 'displacements:' list"),
        ("phonopy:\n  displacements:\n  - atom: 1\n", "no top-level 'displacements:' list"),
        ("displacements: []\n", "no top-level 'displacements:' list"),
        ("displacements:\n- atom: 1\n- displacement: [0.01, 0, 0]\n", ":3: a displacements entry must start"),
        ("displacements:\n- atom: 0\n", ":2: a displacements entry must start"),
        ("displacements:\n- atom: 1\ndisplacements:\n- atom: 2\n", ":3: a second top-level 'displacements:' list"),
        ("displacements:\n  atom: 1\n", ":2: 'displacements:' is not a block list"),
    ],
    ids=["absent", "nested", "flow", "entry", "atom-zero", "twice", "mapping"],
)
def test_an_unsupported_file_is_refused(tmp_path: Path, text: str, message: str) -> None:
    with pytest.raises(ValueError, match=re.escape(message)) as raised:
        read_displacement_count(_write(tmp_path, text))
    # Reasons name the file, not where it is, like the walker's root-relative reasons.
    assert str(raised.value).startswith("phonopy_disp.yaml") and str(tmp_path) not in str(raised.value)


def test_displacement_directories_are_found_in_number_order(tmp_path: Path) -> None:
    for name in ("disp-010", "disp-2", "disp-001", "disp-x", "Disp-003", "disp-004.bak", "other"):
        (tmp_path / name).mkdir()
    (tmp_path / "disp-005").write_text("a file, not a directory\n", encoding="utf-8")
    assert displacement_directories(tmp_path) == ("disp-001", "disp-2", "disp-010")
    assert displacement_directories(tmp_path / "other") == ()
