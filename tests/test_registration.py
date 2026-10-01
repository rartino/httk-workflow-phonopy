"""phonopy registers its collector, and no code, through the ``codes`` registry tier, with its citation."""

import subprocess
import sys

import httk.core  # noqa: F401  (importing httk.core runs registry discovery)
from httk.core.register import collector_support, known_codes, known_collectors


def test_the_collector_is_registered_without_a_code() -> None:
    assert "phonopy.calculation.vasp" in known_collectors()
    assert "phonopy" not in known_codes()
    path = collector_support("phonopy.calculation.vasp").path()
    assert path.name == "phonopy-vasp" and (path / "httk_workflow.toml").is_file()


def test_the_phonopy_credit_is_registered_on_import() -> None:
    script = """
from httk.core import credits
assert "Phonon calculations with phonopy" not in credits.entries()
import httk.codes.phonopy
references = credits.entries()["Phonon calculations with phonopy"]
assert [reference.doi for reference in references] == [
    "10.1088/1361-648X/acd831",
    "10.1016/j.scriptamat.2015.07.021",
]
"""
    subprocess.run([sys.executable, "-c", script], check=True)
