"""Shared test configuration and a synthetic VASP+phonopy finite-displacement tree.

The tree is a CsCl unit cell (``POSCAR``) with phonopy's ``phonopy_disp.yaml``
for a 2x2x2 supercell, whose two symmetry-distinct displacements (one Cs, one
Cl atom) are VASP static runs in ``disp-001`` and ``disp-002``.
"""

from pathlib import Path

import pytest


@pytest.fixture(autouse=True)
def _isolated_httk_config(tmp_path_factory: pytest.TempPathFactory, monkeypatch: pytest.MonkeyPatch) -> None:
    """Give every test its own httk config and data home, so no workspace registry leaks between tests."""

    monkeypatch.setenv("HTTK_CONFIG_HOME", str(tmp_path_factory.mktemp("httk-config")))
    monkeypatch.setenv("HTTK_DATA_HOME", str(tmp_path_factory.mktemp("httk-store")))


A = 4.12  # the CsCl lattice parameter, in Å
SPECIES = (("Cs", 132.905452, (0.0, 0.0, 0.0)), ("Cl", 35.453, (0.5, 0.5, 0.5)))
DISPLACED_ATOMS = (1, 9)  # the first Cs and the first Cl of the supercell, 1-based
STEP = 0.01  # phonopy's default displacement amplitude, in Å
ENERGIES = (-43.21, -43.20)

# A VASP OUTCAR excerpt of a static run, with the lines the OUTCAR reader parses.
_OUTCAR = (
    " vasp.6.4.1 18Apr23 (build Jan 01 2024 00:00:00) complex\n"
    "   NELM   =     60;   NELMIN=  2; NELMDL= -5\n"
    "   NSW    =      0    number of steps for IOM\n"
    "   EDIFF  =   0.100E-05   stopping-criterion for ELM\n"
    "   EDIFFG = -.200E-01   stopping-criterion for IOM\n"
    "   IBRION =     -1    ionic relax: 0-MD;1-quasi-New;2-CG\n"
    "   maximum number of plane-waves:    1234\n"
    "   FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)\n"
    "   free  energy   TOTEN  =       {energy:.8f} eV\n"
    "   energy  without entropy=      {energy:.8f}  energy(sigma->0) =      {energy:.8f}\n"
    " General timing and accounting informations for this job:\n"
)


def _vector(values: tuple[float, ...], digits: int = 15) -> str:
    return "[ " + ", ".join(f"{value:{digits + 4}.{digits}f}" for value in values) + " ]"


def _supercell() -> list[tuple[str, float, tuple[float, float, float], int]]:
    """Return phonopy's 2x2x2 supercell: (symbol, mass, fractional position, unit-cell atom)."""

    atoms = []
    for index, (symbol, mass, (x, y, z)) in enumerate(SPECIES, start=1):
        for k in range(2):
            for j in range(2):
                for i in range(2):
                    atoms.append((symbol, mass, ((x + i) / 2, (y + j) / 2, (z + k) / 2), index))
    return atoms


def _cell_lines(name: str, length: float, points: list[tuple[str, float, tuple[float, float, float], int]]) -> str:
    lines = [f"{name}:", "  lattice:"]
    for axis, label in enumerate("abc"):
        row = tuple(length if column == axis else 0.0 for column in range(3))
        lines.append(f"  - {_vector(row)} # {label}")
    lines.append("  points:")
    for number, (symbol, mass, position, reduced_to) in enumerate(points, start=1):
        lines.append(f"  - symbol: {symbol} # {number}")
        lines.append(f"    coordinates: {_vector(position)}")
        lines.append(f"    mass: {mass:.6f}")
        if name != "primitive_cell":
            lines.append(f"    reduced_to: {reduced_to}")
    if name == "primitive_cell":
        lines.append("  reciprocal_lattice: # without 2pi")
        for axis, label in enumerate("abc"):
            row = tuple(1 / length if column == axis else 0.0 for column in range(3))
            lines.append(f"  - {_vector(row)} # {label}*")
    return "\n".join(lines) + "\n"


def phonopy_disp_yaml(
    displaced: tuple[int, ...] = DISPLACED_ATOMS, *, forces: bool = False, random: bool = False
) -> str:
    """Return a ``phonopy_disp.yaml`` laid out as phonopy 2 writes it.

    :param displaced: The 1-based supercell atoms displaced, one entry each.
    :param forces: Add a ``forces`` list to each entry, as ``phonopy_params.yaml`` has.
    :param random: Write a random-displacement dataset (``phonopy -d --rd N``)
        instead, one ``- # N`` entry per element of *displaced*, each displacing
        every supercell atom.
    """

    unit = [(symbol, mass, position, index) for index, (symbol, mass, position) in enumerate(SPECIES, 1)]
    supercell = [(symbol, mass, position, 1 + 8 * (index - 1)) for symbol, mass, position, index in _supercell()]
    text = (
        "phonopy:\n"
        '  version: "2.26.6"\n'
        "  frequency_unit_conversion_factor: 15.633302\n"
        "  symmetry_tolerance: 1.00000e-05\n"
        "  configuration:\n"
        '    cell_filename: "POSCAR"\n'
        '    create_displacements: ".true."\n'
        '    dim: "2 2 2"\n'
        "  calculator: vasp\n"
        "  physical_unit:\n"
        '    atomic_mass: "AMU"\n'
        '    length: "angstrom"\n'
        '    force_constants: "eV/angstrom^2"\n'
        "\n"
        "space_group:\n"
        '  type: "Pm-3m"\n'
        "  number: 221\n"
        '  Hall_symbol: "-P 4 2 3"\n'
        "\n"
        "primitive_matrix:\n"
        + "".join(f"- {_vector(tuple(1.0 if c == r else 0.0 for c in range(3)))}\n" for r in range(3))
        + "\n"
        "supercell_matrix:\n"
        + "".join("- [ " + ", ".join(f"{2 if c == r else 0:2d}" for c in range(3)) + " ]\n" for r in range(3))
        + "\n"
        + _cell_lines("primitive_cell", A, [(s, m, p, 1) for s, m, p, _ in unit])
        + "\n"
        + _cell_lines("unit_cell", A, unit)
        + "\n"
        + _cell_lines("supercell", 2 * A, supercell)
        + "\n"
        "displacements:\n"
    )
    for number, atom in enumerate(displaced, start=1):
        if random:
            text += f"- # {number:4d}\n" + "".join(
                f"  - displacement: # {n}\n      {_vector((0.0, STEP, 0.0), 16)}\n" for n in range(1, 17)
            )
            continue
        text += f"- atom: {atom}\n  displacement:\n    {_vector((STEP, 0.0, 0.0), 16)}\n"
        if forces:
            text += "  forces:\n" + "".join(f"  - {_vector((0.0, 0.0, 0.0), 10)} # {n}\n" for n in range(1, 17))
    return text


def poscar(displaced: int | None = None) -> str:
    """Return the unit-cell POSCAR, or the supercell POSCAR with one atom displaced (1-based; 0 displaces none)."""

    if displaced is None:
        rows = [position for _symbol, _mass, position in SPECIES]
        length, counts, header = A, "1 1", "CsCl"
    else:
        rows = [position for _symbol, _mass, position, _index in _supercell()]
        if displaced:
            x, y, z = rows[displaced - 1]
            rows[displaced - 1] = (x + STEP / (2 * A), y, z)
        length, counts, header = 2 * A, "8 8", "generated by phonopy"
    lattice = "".join(" ".join(f"{length if c == r else 0.0:21.15f}" for c in range(3)) + "\n" for r in range(3))
    sites = "".join(" ".join(f"{value:19.16f}" for value in row) + "\n" for row in rows)
    return f"{header}\n   1.0\n{lattice}Cs Cl\n{counts}\nDirect\n{sites}"


def vasp_static(directory: Path, displaced: int, energy: float | None, *, encut: int = 520) -> Path:
    """Write a finished VASP static run of one displaced supercell; *energy* ``None`` leaves it unfinished."""

    directory.mkdir(parents=True, exist_ok=True)
    outcar = _OUTCAR.format(energy=0.0 if energy is None else energy)
    if energy is None:
        outcar = outcar.split("   FREE ENERGIE")[0]
    (directory / "INCAR").write_text(
        f"PREC = Accurate\nENCUT = {encut}\nIBRION = -1\nNSW = 0\nEDIFF = 1E-6\n", encoding="utf-8"
    )
    (directory / "POSCAR").write_text(poscar(displaced), encoding="utf-8")
    (directory / "OUTCAR").write_text(outcar, encoding="utf-8")
    return directory


def phonopy_tree(directory: Path, *, energies: tuple[float | None, ...] = ENERGIES, encut: int = 520) -> Path:
    """Write the phonopy calculation into *directory*, one displacement run per entry of *energies*."""

    directory.mkdir(parents=True, exist_ok=True)
    (directory / "phonopy_disp.yaml").write_text(phonopy_disp_yaml(), encoding="utf-8")
    (directory / "POSCAR").write_text(poscar(), encoding="utf-8")
    (directory / "SPOSCAR").write_text(poscar(0), encoding="utf-8")
    for number, (atom, energy) in enumerate(zip(DISPLACED_ATOMS, energies, strict=False), start=1):
        (directory / f"POSCAR-{number:03d}").write_text(poscar(atom), encoding="utf-8")
        vasp_static(directory / f"disp-{number:03d}", atom, energy, encut=encut)
    return directory
