#!/usr/bin/env python3
"""Harmonic phonons of a relaxed structure with phonopy, with the forces from the bound ``forces`` call slot.

This runner computes no forces itself and knows no simulation code. The job's
creator binds the ``forces`` slot (``calls/structure-forces.json``: a structure
in, the forces on its sites out) to a workflow such as ``ase.static`` or
``vasp.static``, and this runner calls it once per displaced supercell:

* ``displace`` standardizes the staged structure to its HPKOT primitive cell with
  ``httk.atomistic.band_path`` (with ``standardize = false`` the given cell is
  kept), makes it right-handed, and orders its sites by species with a stable
  sort; that cell, ``unitcell/POSCAR``, is the cell phonopy diagonalizes, and
  the q-point path of ``unitcell/band_path.json`` refers to its reciprocal
  basis. phonopy builds the diagonal supercell ``n_i = ceil(min_length / w_i)``
  (``w_i`` the cell's width perpendicular to the other two vectors) and its
  displacements (``phonopy_disp.yaml``); each displaced supercell, species
  contiguous because phonopy groups supercell atoms by unit-cell atom, is
  called as ``displacement-<NNN>``.
* ``analyse`` reads each child's recorded forces in displacement order, checks
  their shape, records them in ``phonons.json``, and lets phonopy produce the
  force constants (``FORCE_SETS``), the mesh, the density of states and the band
  structure along the HPKOT path; httk-analyse's ``harmonic_from_phonopy`` gives
  the thermodynamics.

There is no relaxation and no non-analytical correction. Code-specific
consistency across the displacement children (k-points, smearing, remedies)
is the binding's responsibility. ``collect.py`` reads only ``phonons.json``,
``unitcell/POSCAR`` and phonopy's files.
"""

import json
import math
import shutil
from fractions import Fraction
from pathlib import Path
from typing import Any

import httk.atomistic
import httk.core
from httk.workflow import Attempt, Runner

RESULT = "phonons.json"
UNITCELL = Path("unitcell") / "POSCAR"
BAND_PATH = Path("unitcell") / "band_path.json"

run = Runner("phonopy.phonons")


# -- parameters -------------------------------------------------------------------------------------------------------


def real(value: object, what: str, *, positive: bool = True) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{what} must be a number, not {value!r}")
    number = float(value)
    if not math.isfinite(number) or number < 0 or (positive and number == 0):
        raise ValueError(f"{what} must be a finite {'positive' if positive else 'nonnegative'} number, not {value!r}")
    return number


def integer(value: object, what: str, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{what} must be an integer of at least {minimum}, not {value!r}")
    return value


def settings(a: Attempt) -> dict:
    """Return this job's parameters, validated; the manifest declares their types and defaults.

    :raises ValueError: Naming the first invalid parameter.
    """

    temperatures = a.parameter("temperatures")
    if not isinstance(temperatures, list) or not temperatures:
        raise ValueError(f"temperatures must be a nonempty list of numbers, not {temperatures!r}")
    imaginary = a.parameter("imaginary")
    if imaginary not in ("raise", "omit"):
        raise ValueError(f"imaginary must be raise or omit, not {imaginary!r}")
    return {
        "standardize": a.parameter("standardize") is True,
        "min_length": real(a.parameter("min_length"), "min_length"),
        "max_atoms": integer(a.parameter("max_atoms"), "max_atoms", 1),
        "displacement": real(a.parameter("displacement"), "displacement"),
        "mesh_length": real(a.parameter("mesh_length"), "mesh_length"),
        "band_points": integer(a.parameter("band_points"), "band_points", 2),
        "temperatures": [real(t, "a temperature", positive=False) for t in temperatures],
        "cutoff_frequency": real(a.parameter("cutoff_frequency"), "cutoff_frequency", positive=False),
        "imaginary": imaginary,
        "imaginary_tolerance": real(a.parameter("imaginary_tolerance"), "imaginary_tolerance", positive=False),
        "symmetry_tolerance": real(a.parameter("symmetry_tolerance"), "symmetry_tolerance"),
    }


def write_json(path: Path, document: dict) -> None:
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(document, indent=1), encoding="utf-8")
    temporary.replace(path)


def read_result(a: Attempt) -> dict:
    return json.loads((a.workdir / RESULT).read_text(encoding="utf-8"))


def forces_problem(value: object, atoms: int) -> str | None:
    """Return why *value* is not an ``atoms x 3`` array of finite numbers, or ``None`` when it is."""

    if not isinstance(value, list) or len(value) != atoms:
        return (
            f"expected {atoms} rows, one per supercell atom, got {len(value) if isinstance(value, list) else value!r}"
        )
    for index, row in enumerate(value):
        if not isinstance(row, list) or len(row) != 3:
            return f"row {index} is not a 3-vector: {row!r}"
        if any(isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) for x in row):
            return f"row {index} is not finite numbers: {row!r}"
    return None


# -- steps ------------------------------------------------------------------------------------------------------------


@run.step
def displace(a: Attempt) -> None:
    """Fix the phonon unit cell and its q-point path, and call the forces slot once per displaced supercell."""

    try:
        own = settings(a)
    except ValueError as exception:
        a.fail("phonons.input_invalid", str(exception))
        return
    try:
        import numpy
        from phonopy import Phonopy
        from phonopy.interface.vasp import read_vasp, write_vasp
        from phonopy.structure.atoms import PhonopyAtoms
    except ImportError as exception:
        a.fail("phonons.analysis_unavailable", f"phonons need phonopy: {exception} (pip install phonopy)")
        return
    staged = a.payload / "files" / "POSCAR"
    unitcell = a.workdir / UNITCELL
    unitcell.parent.mkdir(exist_ok=True)
    tolerance = own["symmetry_tolerance"]
    try:
        structure = httk.core.load(str(staged), format="vasp-poscar", precision=tolerance)
        result = httk.atomistic.band_path(structure, tolerance=tolerance)
    except (ImportError, ValueError) as exception:
        a.fail("phonons.band_path_failed", f"cannot find the HPKOT band path of the structure: {exception}")
        return
    if own["standardize"]:
        httk.core.save(result.structure, str(unitcell))
    else:
        shutil.copyfile(staged, unitcell)
    given = read_vasp(unitcell)
    lattice = numpy.array(given.cell)
    positions = numpy.array(given.scaled_positions)
    # Right-handed, once on the unit cell so every supercell inherits it: a2 and a3 swap, and with them the
    # fractional coordinate columns, so the crystal stays the same (not its mirror image); the basis change
    # below absorbs the swap.
    if numpy.linalg.det(lattice) < 0:
        lattice, positions = lattice[[0, 2, 1]], positions[:, [0, 2, 1]]
    # Species-contiguous (a stable sort by first appearance): writers such as POSCAR group sites by species, and
    # the forces slot refuses a staged supercell whose site order differs from the one passed.
    symbols = list(given.symbols)
    order = sorted(range(len(symbols)), key=lambda i: symbols.index(symbols[i]))
    cell = PhonopyAtoms(
        symbols=[symbols[i] for i in order],
        cell=lattice,
        scaled_positions=positions[order],
    )
    write_vasp(unitcell, cell)
    lattice_rows: list[list[float]] = cell.cell.tolist()
    path_lattice = [[float(x) for x in row] for row in result.structure.cell.basis.to_floats()]
    # band_path's transform starts from httk's primitive cell, not the unit cell, so the cells are compared: with
    # rows as vectors, L_path = T L_unit, and the unit cell is accepted when T is a unimodular integer matrix.
    # Standardized, T is the identity, or the a2/a3 swap when the handedness was normalized.
    basis_change = basis_change_of(path_lattice, lattice_rows)
    if basis_change is None:
        a.fail(
            "phonons.input_invalid" if not own["standardize"] else "phonons.band_path_failed",
            "the given unit cell is not a basis of its HPKOT lattice, so the q-point path would not refer to it: "
            "pass standardize = true",
            details={"unit_cell": lattice_rows, "hpkot_cell": path_lattice},
        )
        return
    symbol = result.bravais_lattice_extended
    write_json(
        a.workdir / BAND_PATH,
        {
            "convention": "HPKOT",
            "bravais_lattice": result.bravais_lattice,
            "bravais_lattice_extended": symbol,
            "spacegroup": result.spacegroup.it_number,
            "symmetry_tolerance": tolerance,
            # Fractional coordinates of the unit cell's reciprocal basis, transformed from band_path's cell.
            "points": {
                label: [str(v) for v in to_unit_cell(point.to_fractions(), basis_change)]
                for label, point in result.points.items()
            },
            "path": [list(segment) for segment in result.path],
            "transform": [[int(v) for v in row] for row in result.transform.to_fractions()],
            "basis_change": basis_change,
        },
    )

    # The diagonal supercell by perpendicular widths w_i = V / |a_j x a_k|.
    a1, a2, a3 = lattice_rows
    volume = abs(_dot(a1, _cross(a2, a3)))
    widths = [volume / math.sqrt(_dot(c, c)) for c in (_cross(a2, a3), _cross(a3, a1), _cross(a1, a2))]
    # ponytail: diagonal supercells only, as vasp.phonons; atomate2's cubic-like (non-diagonal) search is the upgrade.
    multiples = [math.ceil(own["min_length"] / w - 1e-9) for w in widths]
    atoms = len(symbols) * multiples[0] * multiples[1] * multiples[2]
    if atoms > own["max_atoms"]:
        a.fail(
            "phonons.supercell_too_large",
            f"the {'x'.join(map(str, multiples))} supercell for min_length {own['min_length']} has {atoms} atoms, "
            f"more than max_atoms = {own['max_atoms']}",
            details={"multiples": multiples, "atoms": atoms},
        )
        return
    supercell_matrix = [[n if i == j else 0 for j in range(3)] for i, n in enumerate(multiples)]
    # Any: phonopy types its results Optional until the run_* call that fills them.
    phonon: Any = Phonopy(cell, supercell_matrix=supercell_matrix, primitive_matrix="P")
    phonon.generate_displacements(distance=own["displacement"])
    phonon.save(a.workdir / "phonopy_disp.yaml")
    labels = [f"displacement-{index:03d}" for index in range(len(phonon.supercells_with_displacements))]
    for label, displaced in zip(labels, phonon.supercells_with_displacements, strict=True):
        path = a.workdir / "displacements" / label / "POSCAR"
        path.parent.mkdir(parents=True, exist_ok=True)
        write_vasp(path, displaced)
        try:
            # Read back as an httk structure, so the bound workflow is staged it in its own format.
            supercell = httk.core.load(str(path), format="vasp-poscar", precision=1e-9)
            problem = written_problem(supercell, displaced)
            if problem is not None:
                a.fail("phonons.input_invalid", f"{label}: the written POSCAR is not phonopy's supercell: {problem}")
                return
            a.call("forces", label=label, inputs={"structure": supercell})
        except ValueError as exception:
            a.fail("phonons.call_refused", f"{label}: {exception}", details={"label": label})
            return
    write_json(
        a.workdir / RESULT,
        {
            "workflow": "phonopy.phonons",
            "parameters": own,
            "unit_cell": {
                "symbols": cell.symbols,
                "lattice": lattice_rows,
                "positions": cell.scaled_positions.tolist(),
            },
            "band_path": {"bravais_lattice_extended": symbol, "spacegroup": result.spacegroup.it_number},
            "supercell_matrix": supercell_matrix,
            "supercell_atoms": atoms,
            # Displacement index i is child label displacement-<i>, in phonopy_disp.yaml's order.
            "displacements": [
                {
                    "label": label,
                    "atom": int(entry["number"]),
                    "displacement": [float(x) for x in entry["displacement"]],
                }
                for label, entry in zip(labels, phonon.dataset["first_atoms"], strict=True)
            ],
        },
    )
    a.state.merge({"supercell": "x".join(map(str, multiples)), "displacements": len(labels)})
    a.log.append(
        "note",
        f"standardize = {own['standardize']}: phonons of the {symbol} cell in a {'x'.join(map(str, multiples))} "
        f"supercell of {atoms} atoms; called the forces slot for {len(labels)} displacements",
    )
    # all_terminal: analyse records a failed displacement before failing on it.
    a.gather("analyse", when="all_terminal")


def written_problem(structure: object, displaced: Any) -> str | None:
    """Return how *structure*, read back from phonopy's POSCAR, differs from phonopy's cell, or ``None``.

    The forces come back in the order of the structure passed to the call; this makes that order phonopy's.
    """

    view = httk.atomistic.UnitcellStructureView(structure)  # type: ignore[arg-type]
    by_name = {species.name: tuple(species.chemical_symbols) for species in view.species}
    symbols = [by_name[name] for name in view.species_at_sites]
    if symbols != [(symbol,) for symbol in displaced.symbols]:
        return f"its sites are {symbols}, phonopy's {list(displaced.symbols)}"
    wanted = displaced.scaled_positions.tolist()
    for index, (got, want) in enumerate(zip(view.sites.reduced_coords.to_floats(), wanted, strict=True)):
        if any(abs((x - y) - round(x - y)) > 1e-9 for x, y in zip(got, want, strict=True)):
            return f"site {index} is at {list(got)}, phonopy's at {want}"
    return None


def basis_change_of(path_lattice: list[list[float]], lattice: list[list[float]]) -> list[list[int]] | None:
    """Return the integer ``T`` with ``path_lattice = T lattice`` (rows as vectors), or ``None`` unless it is unimodular."""

    import numpy

    change = numpy.array(path_lattice) @ numpy.linalg.inv(numpy.array(lattice))
    rounded = numpy.rint(change)
    if numpy.abs(change - rounded).max() > 1e-6 or round(abs(numpy.linalg.det(rounded))) != 1:
        return None
    return [[int(x) for x in row] for row in rounded]


def to_unit_cell(point: list[Fraction], basis_change: list[list[int]]) -> list[Fraction]:
    """Return a k-point of the path cell's reciprocal basis in the unit cell's: ``q T^-T`` for ``L_path = T L_unit``.

    Reciprocal rows are ``B = L^-T``, so ``B_path = T^-T B_unit`` and ``k = q B_path = (q T^-T) B_unit``.
    """

    import numpy

    # T is unimodular, so its inverse is an integer matrix too.
    inverse_transposed = numpy.rint(numpy.linalg.inv(numpy.array(basis_change, dtype=float))).astype(int).T
    return [sum((point[i] * int(inverse_transposed[i][j]) for i in range(3)), Fraction(0)) for j in range(3)]


def _dot(u: list[float], v: list[float]) -> float:
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def _cross(u: list[float], v: list[float]) -> list[float]:
    return [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]


@run.step
def analyse(a: Attempt) -> None:
    """Read every displacement child's forces in displacement order and compute the phonons from them.

    Everything read is persisted (the children's ``outputs.json``, this workdir), so a retried analysis reads the
    same forces and gives the same result.
    """

    own = settings(a)
    result = read_result(a)
    atoms = result["supercell_atoms"]
    records: list[dict] = []
    unfinished: dict[str, dict] = {}
    for entry in result["displacements"]:
        label = entry["label"]
        child = a.children.get(label)
        if child is None:
            unfinished[label] = {"label": label, "job_id": None, "kind": "missing", "failure": None}
        elif not child.succeeded:
            failure = None if child.failure is None else {"code": child.failure.code, "message": child.failure.message}
            unfinished[label] = {"label": label, "job_id": child.job_id, "kind": child.kind, "failure": failure}
    if unfinished:
        result["children"] = list(unfinished.values())
        write_json(a.workdir / RESULT, result)
        failed = [label for label, entry in unfinished.items() if entry["kind"] != "missing"]
        a.fail(
            "phonons.child_failed" if failed else "phonons.child_missing",
            f"displacements {', '.join(failed or unfinished)} "
            + ("did not succeed" if failed else "have no child job")
            + "; phonons.json records them",
            details={"unfinished": result["children"]},
        )
        return
    for entry in result["displacements"]:
        label = entry["label"]
        child = a.children[label]
        try:
            forces = a.child_outputs(child, call="forces")["forces"]
        except ValueError as exception:
            a.fail("phonons.forces_unavailable", f"{label}: {exception}", details={"label": label})
            return
        problem = forces_problem(forces, atoms)
        if problem is not None:
            a.fail("phonons.forces_invalid", f"the forces of {label} are invalid: {problem}", details={"label": label})
            return
        # ponytail: no remedy histories; their location is code-specific (VASP's), so a binding records its own.
        records.append({"label": label, "job_key": child.job_key, "job_id": child.job_id, "forces": forces})
    result["children"] = records
    try:
        import numpy
        import phonopy
        from httk.analyse.integrations.phonopy import harmonic_from_phonopy
        from phonopy.file_IO import write_FORCE_SETS
    except ImportError as exception:
        a.fail(
            "phonons.analysis_unavailable",
            f"the phonon analysis needs phonopy and httk-analyse: {exception} (pip install phonopy httk-analyse)",
        )
        return
    phonon: Any = phonopy.load(a.workdir / "phonopy_disp.yaml", produce_fc=False, is_nac=False)
    if len(phonon.dataset["first_atoms"]) != len(records):
        a.fail("phonons.inconsistent_children", "phonopy_disp.yaml does not list the displacements that were called")
        return
    phonon.forces = numpy.array([record["forces"] for record in records])
    phonon.produce_force_constants()
    write_FORCE_SETS(phonon.dataset, a.workdir / "FORCE_SETS")
    phonon.run_mesh(own["mesh_length"])
    frequencies = numpy.asarray(phonon.mesh.frequencies)
    imaginary = frequencies < -own["imaginary_tolerance"]
    result["mesh"] = {
        "length": own["mesh_length"],
        "numbers": [int(n) for n in phonon.mesh.mesh_numbers],
        "irreducible_qpoints": len(frequencies),
    }
    result["imaginary_modes"] = {
        "tolerance": own["imaginary_tolerance"],
        "count": int(imaginary.sum()),
        "lowest_frequency": float(frequencies.min()),
    }
    try:
        thermodynamics = harmonic_from_phonopy(
            phonon,
            own["temperatures"],
            zero_modes="omit",
            imaginary=own["imaginary"],
            cutoff_frequency=own["cutoff_frequency"],
        )
    except ValueError as exception:
        result["thermodynamics"] = None
        result["thermodynamics_omitted"] = str(exception)
        a.log.append("note", f"no thermodynamics (imaginary = {own['imaginary']}): {exception}")
    else:
        result["thermodynamics"] = {name: getattr(thermodynamics, name) for name in thermodynamics.__dataclass_fields__}
    if imaginary.any():
        a.log.append(
            "note",
            f"{int(imaginary.sum())} imaginary modes on the mesh below -{own['imaginary_tolerance']} THz, "
            f"the lowest {frequencies.min():.4f} THz",
        )
    # phonopy sizes the DOS grid by the spectrum's width, which a dispersionless (Einstein) spectrum lacks: a
    # 0.01 THz Gaussian then gives it one.
    phonon.run_total_dos(sigma=None if numpy.ptp(frequencies) > 1e-6 else 0.01)
    result["dos"] = {
        "frequencies": [float(x) for x in phonon.total_dos.frequency_points],
        "density": [float(x) for x in phonon.total_dos.dos],
    }
    path = json.loads((a.workdir / BAND_PATH).read_text(encoding="utf-8"))
    points = {label: [float(Fraction(x)) for x in point] for label, point in path["points"].items()}
    count = own["band_points"]
    segments = [
        [[s + (e - s) * i / (count - 1) for s, e in zip(points[start], points[end], strict=True)] for i in range(count)]
        for start, end in path["path"]
    ]
    phonon.run_band_structure(segments)
    band = phonon.band_structure
    stacked = numpy.concatenate([numpy.asarray(f) for f in band.frequencies])
    result["band_structure"] = {
        "qpoints": [list(map(float, q)) for segment in band.qpoints for q in segment],
        "qpoint_label_indices": [i * count + j for i in range(len(segments)) for j in (0, count - 1)],
        "qpoint_labels": [label for segment in path["path"] for label in segment],
        # Sorted at each q-point, as atomate2.phonons' collect does: frequencies[mode][q].
        "frequencies": numpy.sort(stacked, axis=1).T.tolist(),
    }
    write_json(a.workdir / RESULT, result)
    headline = {
        "imaginary_modes": result["imaginary_modes"]["count"],
        "lowest_frequency": result["imaginary_modes"]["lowest_frequency"],
        "zero_point_energy": None
        if result["thermodynamics"] is None
        else result["thermodynamics"]["zero_point_energy"],
    }
    a.state.merge(headline)
    a.log.append("note", f"phonons: {headline}")
    a.succeed()


if __name__ == "__main__":
    raise SystemExit(run.main())
