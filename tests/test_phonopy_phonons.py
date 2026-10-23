"""The ``phonopy.phonons`` workflow end to end, with an ASE EMT callee bound to its ``forces`` call slot.

The callee is an in-test copy of the ``ase.static`` example of *workflows-examples*,
so this repository's tests need no sibling checkout; when that checkout exists the
fcc Cu case also runs against it. Each run uses a real local manager, and each
result is compared with phonopy driven directly with ASE's EMT on the same cell.
"""

import json
import shutil
import tomllib
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest

pytest.importorskip("phonopy")
pytest.importorskip("ase")
pytest.importorskip("httk.atomistic")
pytest.importorskip("httk.analyse")

REPO_ROOT = Path(__file__).resolve().parents[1]
PACKAGE = REPO_ROOT / "workflows" / "phonopy-phonons"
ASE_STATIC = REPO_ROOT.parent / "workflows-examples" / "ase-static"
STRUCTURE = "https://schemas.optimade.org/defs/v1.3/entrytypes/optimade/structures"
FORCES = "https://schemas.httk.org/defs/v0.1/properties/core/cartesian_site_forces"

CU = """fcc Cu
1.0
0.0 1.795 1.795
1.795 0.0 1.795
1.795 1.795 0.0
Cu
1
Direct
0.0 0.0 0.0
"""
# L1_2 Cu3Au with its sites interleaved by species: Cu, Au, Cu, Cu.
CU3AU_INTERLEAVED = """Cu3Au
1.0
3.75 0.0 0.0
0.0 3.75 0.0
0.0 0.0 3.75
Cu Au Cu Cu
1 1 1 1
Direct
0.5 0.5 0.0
0.0 0.0 0.0
0.5 0.0 0.5
0.0 0.5 0.5
"""
# Small supercells (2x2x2 of both cells) and a coarse mesh keep the runs to seconds.
PARAMETERS = {"min_length": 4.0, "mesh_length": 10.0, "band_points": 5, "temperatures": [0.0, 300.0]}

_CALLEE_RUN = '''#!/usr/bin/env python3
import importlib
import json
import os

from httk.workflow import Runner

run = Runner("tests.ase.static")


@run.step
def evaluate(a):
    import ase.io

    marker = a.parameter("fail_once")
    if marker and os.path.exists(marker):
        os.remove(marker)
        a.fail("tests.forced", "failed once on purpose")
        return
    module, _, name = a.parameter("calculator").partition(":")
    atoms = ase.io.read(a.payload / "files" / "POSCAR", format="vasp")
    atoms.calc = getattr(importlib.import_module(module), name)()
    results = {"total_energy": float(atoms.get_potential_energy()), "forces": atoms.get_forces().tolist()}
    (a.workdir / "results.json").write_text(json.dumps(results), encoding="utf-8")
    a.succeed()


if __name__ == "__main__":
    raise SystemExit(run.main())
'''

_CALLEE_COLLECT = """import json

from httk.core.data_records import TotalEnergyRecord
from httk.core.property_records import CartesianSiteForcesRecord


def collect(record):
    results = json.loads(record.result_file("results.json").read_text(encoding="utf-8"))
    return {
        "total_energy": TotalEnergyRecord(total_energy=results["total_energy"]),
        "cartesian_site_forces": CartesianSiteForcesRecord.from_value(results["forces"]),
    }
"""


def _callee(root: Path, *, forces_of: str = "initial_structure", procs: int | None = None) -> Path:
    """Write the ASE single-point callee package; its forces are a product of *forces_of*."""

    manifest = f"""[workflow]
name = "tests.ase.static"

[workflow.runner]
entry = "run.py"
steps = ["evaluate"]

[workflow.inputs.structure]
destination = "POSCAR"
role = "initial_structure"
entry_type = "structures"
ref = "{STRUCTURE}"
required = true

[workflow.parameters.calculator]
type = "string"
default = "ase.calculators.emt:EMT"

[workflow.parameters.fail_once]
type = "string"
default = ""

[workflow.outputs.total_energy]
entry_type = "records"
role = "total_energy"
ref = "https://schemas.httk.org/defs/v0.1/properties/core/total_energy"
product_of = "initial_structure"

[workflow.outputs.cartesian_site_forces]
entry_type = "records"
role = "cartesian_site_forces"
ref = "{FORCES}"
product_of = "{forces_of}"

[workflow.collect]
file = "collect.py"
"""
    if forces_of != "initial_structure":
        manifest += f'\n[workflow.outputs.{forces_of}]\nentry_type = "structures"\nproduct_of = "initial_structure"\n'
    if procs is not None:
        manifest += f"\n[workflow.resources]\nprocs = {procs}\n"
    root.mkdir(parents=True)
    (root / "httk_workflow.toml").write_text(manifest, encoding="utf-8")
    (root / "collect.py").write_text(_CALLEE_COLLECT, encoding="utf-8")
    (root / "run.py").write_text(_CALLEE_RUN, encoding="utf-8")
    (root / "run.py").chmod(0o755)
    return root


@pytest.fixture
def installed_plugin(tmp_path_factory: pytest.TempPathFactory) -> Iterator[None]:
    """Install this repository's ``httk_plugin.toml`` plugin, as ``httk plugin install`` does."""

    from httk.core.plugins import install_plugin
    from httk.workflow.packages import _reset_plugin_workflow_cache

    source = tmp_path_factory.mktemp("plugin-source") / "httk-workflow-phonopy"
    shutil.copytree(REPO_ROOT / "workflows", source / "workflows", ignore=shutil.ignore_patterns("__pycache__"))
    shutil.copy2(REPO_ROOT / "httk_plugin.toml", source)
    install_plugin(source)
    _reset_plugin_workflow_cache()
    yield
    _reset_plugin_workflow_cache()


def _workspace(tmp_path: Path):
    from httk.workflow import Workspace

    return Workspace.initialize(tmp_path / "workspace", durable=False)


def _install(workspace, source: Path) -> str:
    """Install the workflow package at *source* into *workspace* and return its id: a slot binds no paths."""

    from httk.workflow import _kernel, _store

    owner = _kernel.register_owner(workspace, kind="cli", label="test", allocation=None, advertised={})
    try:
        return _store.install(workspace, owner, source).id
    finally:
        owner.close()


def _new_job(
    workspace,
    tmp_path: Path,
    poscar: str,
    callee: Path | str,
    workflow: object = PACKAGE,
    call_parameters: dict | None = None,
    **parameters,
):
    """Submit one job binding *callee* (a package directory, installed first, or an installed id) to ``forces``."""

    from httk.workflow.scaffold import new_job

    (tmp_path / "POSCAR").write_text(poscar, encoding="utf-8")
    return new_job(
        workspace,
        workflow,  # type: ignore[arg-type]
        inputs={"structure": tmp_path / "POSCAR"},
        parameters=PARAMETERS | parameters,
        calls={"forces": _install(workspace, callee) if isinstance(callee, Path) else callee},
        call_parameters={"forces": call_parameters or {}},
        install=True,
    )


def _run(workspace, **manager) -> None:
    from httk.workflow import TaskManager

    with TaskManager(workspace, heartbeat_interval=0.01, **manager) as task_manager:
        task_manager.run_until_idle(timeout=300.0)


def _record(workspace, job_id: str):
    from httk.workflow.collecting import job_records

    [record] = [r for r in job_records(workspace, states=("succeeded", "failed")) if r.job_id == job_id]
    return record


def _succeeded(workspace, job_id: str):
    record = _record(workspace, job_id)
    assert record.state == "succeeded", record.failure
    return record


def _direct(record) -> tuple[dict, list[list[float]], list]:
    """Return the job's ``phonons.json``, phonopy's frequencies on its q-points and the EMT forces, computed here."""

    import numpy
    from ase import Atoms
    from ase.calculators.emt import EMT
    from phonopy import Phonopy
    from phonopy.interface.vasp import read_vasp

    result = json.loads(record.result_file("phonons.json").read_text(encoding="utf-8"))
    phonon: Any = Phonopy(
        read_vasp(record.result_file("unitcell/POSCAR")),
        supercell_matrix=result["supercell_matrix"],
        primitive_matrix="P",
    )
    phonon.generate_displacements(distance=result["parameters"]["displacement"])
    forces = []
    for cell in phonon.supercells_with_displacements:
        atoms = Atoms(cell.symbols, cell=cell.cell, scaled_positions=cell.scaled_positions, pbc=True)
        atoms.calc = EMT()
        forces.append(atoms.get_forces())
    phonon.forces = numpy.array(forces)
    phonon.produce_force_constants()
    phonon.run_qpoints(result["band_structure"]["qpoints"])
    expected = numpy.sort(phonon.qpoints.frequencies, axis=1).T.tolist()
    return result, expected, forces


def _assert_matches_direct(record) -> dict:
    import numpy

    result, expected, forces = _direct(record)
    # Row by row: any permutation of the sites between phonopy, the staged file and the forces shows here.
    assert len(result["children"]) == len(forces)
    for child, direct in zip(result["children"], forces, strict=True):
        assert numpy.allclose(child["forces"], direct, atol=1e-8, rtol=0), child["label"]
    # Compared as signed squares (eigenvalues): a zero mode's frequency is the square root of rounding noise.
    got, want = numpy.array(result["band_structure"]["frequencies"]), numpy.array(expected)
    assert numpy.allclose(numpy.sign(got) * got**2, numpy.sign(want) * want**2, atol=1e-6, rtol=0)
    return result


def _collected(workspace, job_id: str) -> dict:
    from httk.workflow.collecting import collect

    [collected] = [c for c in collect(workspace, allow_job_collector=True, fail_fast=True) if c.record.job_id == job_id]
    return dict(collected.outputs)


def _fcc_cu(tmp_path: Path, callee: Path | str, workflow: object = PACKAGE) -> None:
    workspace = _workspace(tmp_path)
    job = _new_job(workspace, tmp_path, CU, callee, workflow)
    _run(workspace)
    record = _succeeded(workspace, job.job_id)
    result = _assert_matches_direct(record)
    gamma = result["band_structure"]["qpoints"].index([0.0, 0.0, 0.0])
    acoustic = [mode[gamma] for mode in result["band_structure"]["frequencies"]][:3]
    assert all(abs(f) < 0.05 for f in acoustic), acoustic
    assert len(result["children"]) == len(result["displacements"]) == 1
    outputs = _collected(workspace, job.job_id)
    assert {
        "phonon_cell",
        "phonon_band_structure",
        "phonon_density_of_states",
        "vibrational_thermodynamics",
        "zero_point_energy",
        "phonons_result",
        "phonopy_disp",
        "force_sets",
    } <= set(outputs)


def test_fcc_cu_phonons_from_the_emt_forces_slot(tmp_path: Path, installed_plugin: None) -> None:
    _fcc_cu(tmp_path, _callee(tmp_path / "callee"), "phonopy.phonons")


@pytest.mark.skipif(not ASE_STATIC.is_dir(), reason="needs the workflows-examples checkout beside this repository")
def test_fcc_cu_phonons_from_the_ase_static_example(tmp_path: Path) -> None:
    _fcc_cu(tmp_path, ASE_STATIC)


@pytest.mark.parametrize("standardize", [True, False], ids=["standardized", "given-cell"])
def test_interleaved_species_are_sorted_before_phonopy(tmp_path: Path, standardize: bool) -> None:
    workspace = _workspace(tmp_path)
    job = _new_job(workspace, tmp_path, CU3AU_INTERLEAVED, _callee(tmp_path / "callee"), standardize=standardize)
    _run(workspace)
    record = _succeeded(workspace, job.job_id)
    result = _assert_matches_direct(record)
    assert result["unit_cell"]["symbols"] == ["Cu", "Cu", "Cu", "Au"]
    assert result["supercell_atoms"] == 32 and len(result["children"]) == len(result["displacements"])


def test_a_retried_analysis_rereads_the_same_forces(tmp_path: Path) -> None:
    from httk.core.cli import CLIContext
    from httk.core.identity import ensure_identity_key, write_identity_config
    from httk.workflow.collecting import job_records
    from httk.workflow.registry import register_workspace
    from httk.workflow.workflow_cli import command

    workspace = _workspace(tmp_path)
    # The callee asks for two processors: a one-processor manager runs only displace.
    job = _new_job(workspace, tmp_path, CU, _callee(tmp_path / "callee", procs=2))
    _run(workspace, resources={"procs": 1})
    assert not list(job_records(workspace, states=("succeeded", "failed")))
    # The children alone (their placement subtree), then one child's recorded forces disappear for a while.
    children_placement = f"{job.placement}/{job.job_key}.children"
    _run(workspace, resources={"procs": 2}, placement_prefixes=(children_placement,))
    [child] = list(job_records(workspace, states=("succeeded",)))
    outputs = child.payload / ".httk-job" / "outputs.json"
    recorded = json.loads(outputs.read_text(encoding="utf-8"))["outputs"]["forces"]["value"]
    outputs.rename(outputs.with_suffix(".away"))
    _run(workspace, resources={"procs": 2})
    failed = _record(workspace, job.job_id)
    assert failed.state == "failed" and failed.failure is not None
    assert failed.failure.code == "phonons.forces_unavailable" and "displacement-000" in failed.failure.message

    outputs.with_suffix(".away").rename(outputs)
    write_identity_config({"identities": {"op": {"name": "Op", "email": "op@example.test"}}, "default_identity": "op"})
    ensure_identity_key("op")
    register_workspace("phonons", str(workspace.root))
    arguments = ["job", "request", "continue", "--workspace", "phonons", "--reason", "outputs back", job.job_id]
    assert command(arguments, CLIContext("httk", tmp_path)) == 0
    _run(workspace, resources={"procs": 2})
    result = _assert_matches_direct(_succeeded(workspace, job.job_id))
    assert [entry["forces"] for entry in result["children"]] == [recorded]


def test_a_forces_workflow_of_another_structure_does_not_fit_the_slot(tmp_path: Path) -> None:
    workspace = _workspace(tmp_path)
    callee = _callee(tmp_path / "callee", forces_of="relaxed_structure")
    with pytest.raises(ValueError, match="call slot 'forces' bound to .*interface output 'forces'"):
        _new_job(workspace, tmp_path, CU, callee)


VASP_STATIC = REPO_ROOT.parent / "workflows-vasp" / "vasp-static"
#: The forces slot's ``default_parameters``, as the manifest declares them.
DEFAULT_PARAMETERS = tomllib.loads((PACKAGE / "httk_workflow.toml").read_text(encoding="utf-8"))["workflow"]["calls"][
    "forces"
]["default_parameters"]


def _vasp_static_stub(root: Path) -> Path:
    """Write a package named ``vasp.static`` declaring what the real one does (and no parameters); it never runs."""

    root.mkdir(parents=True)
    (root / "httk_workflow.toml").write_text(
        f"""[workflow]
name = "vasp.static"

[workflow.runner]
entry = "run.py"
steps = ["prepare"]

[workflow.inputs.structure]
destination = "POSCAR"
role = "initial_structure"
entry_type = "structures"
ref = "{STRUCTURE}"
required = true

[workflow.outputs.total_energy]
entry_type = "records"
role = "total_energy"
ref = "https://schemas.httk.org/defs/v0.1/properties/core/total_energy"
product_of = "initial_structure"

[workflow.outputs.cartesian_site_forces]
entry_type = "records"
role = "cartesian_site_forces"
ref = "{FORCES}"
product_of = "initial_structure"
""",
        encoding="utf-8",
    )
    (root / "run.py").write_text("#!/usr/bin/env python3\nraise SystemExit(1)\n", encoding="utf-8")
    (root / "run.py").chmod(0o755)
    return root


def _default_binding(tmp_path: Path, vasp_static: Path | None, **call_parameters) -> tuple[str | None, dict]:
    """Create a job that binds nothing, with *vasp_static* installed first when given.

    :return: The installed id of *vasp_static* (or ``None``) and the job's recorded ``forces`` binding.
    """

    from httk.workflow.scaffold import new_job

    workspace = _workspace(tmp_path)
    installed = None if vasp_static is None else _install(workspace, vasp_static)
    _install(workspace, PACKAGE)
    (tmp_path / "POSCAR").write_text(CU, encoding="utf-8")
    job = new_job(
        workspace,
        "phonopy.phonons",
        inputs={"structure": tmp_path / "POSCAR"},
        call_parameters={"forces": call_parameters} if call_parameters else None,
    )
    return installed, json.loads((job.payload / "job.json").read_text(encoding="utf-8"))["calls"]["forces"]


def test_an_unbound_slot_takes_vasp_static_with_its_phonon_safe_parameters(tmp_path: Path) -> None:
    installed, binding = _default_binding(tmp_path, _vasp_static_stub(tmp_path / "vasp-static"))
    assert binding == {"workflow": installed, "parameters": DEFAULT_PARAMETERS}
    assert DEFAULT_PARAMETERS["incar_tags"]["EDIFF"] == "1E-7" and DEFAULT_PARAMETERS["kpoint_density"] == 20


def test_a_call_parameter_overrides_one_default_parameter(tmp_path: Path) -> None:
    _installed, binding = _default_binding(tmp_path, _vasp_static_stub(tmp_path / "vasp-static"), kpoint_density=30)
    assert binding["parameters"] == DEFAULT_PARAMETERS | {"kpoint_density": 30}


def test_another_bound_workflow_gets_no_default_parameters(tmp_path: Path) -> None:
    workspace = _workspace(tmp_path)
    job = _new_job(workspace, tmp_path, CU, _callee(tmp_path / "callee"))
    binding = json.loads((job.payload / "job.json").read_text(encoding="utf-8"))["calls"]["forces"]
    assert binding["parameters"] == {}


def test_the_default_must_be_installed(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="call slot 'forces' of .* defaults to vasp.static, which is not installed"):
        _default_binding(tmp_path, None)


@pytest.mark.skipif(not VASP_STATIC.is_dir(), reason="needs the workflows-vasp checkout beside this repository")
def test_the_real_vasp_static_fits_the_slot_with_the_default_parameters(tmp_path: Path) -> None:
    pytest.importorskip("httk.codes.vasp")
    installed, binding = _default_binding(tmp_path, VASP_STATIC)
    assert binding == {"workflow": installed, "parameters": DEFAULT_PARAMETERS}


def _shifted_cu3au(handedness: str) -> str:
    """Return L1_2 Cu3Au shifted off its inversion centre, in a right- or left-handed (a2, a3 swapped) cell."""

    sites = [(0.6, 0.7, 0.3), (0.6, 0.2, 0.8), (0.1, 0.7, 0.8), (0.1, 0.2, 0.3)]
    rows = ["3.75 0.0 0.0", "0.0 3.75 0.0", "0.0 0.0 3.75"]
    if handedness == "left":
        rows, sites = [rows[0], rows[2], rows[1]], [(x, z, y) for x, y, z in sites]
    lines = ["Cu3Au", "1.0", *rows, "Cu Au", "3 1", "Direct", *(" ".join(map(str, site)) for site in sites)]
    return "\n".join(lines) + "\n"


def _cartesian_sites(poscar: str) -> list[tuple[str, tuple[float, ...]]]:
    """Return the sites of a POSCAR as (symbol, fractional position mod 1 in the cubic 3.75 Angstrom frame)."""

    import numpy

    lines = poscar.splitlines()
    lattice = numpy.array([[float(x) for x in line.split()] for line in lines[2:5]])
    symbols = [s for s, n in zip(lines[5].split(), lines[6].split(), strict=True) for _ in range(int(n))]
    fractional = numpy.array([[float(x) for x in line.split()[:3]] for line in lines[8 : 8 + len(symbols)]])
    # Cartesian, then fractional in the axis frame: equal sets mean the same crystal up to lattice translations.
    framed = (fractional @ lattice / 3.75) % 1.0
    return sorted((s, tuple(round(x, 6) % 1.0 for x in row)) for s, row in zip(symbols, framed, strict=True))


def test_a_left_handed_cell_keeps_its_crystal_and_its_phonons(tmp_path: Path) -> None:
    import numpy

    workspace = _workspace(tmp_path)
    callee = _install(workspace, _callee(tmp_path / "callee"))
    jobs = {}
    for handedness in ("right", "left"):
        (tmp_path / handedness).mkdir()
        jobs[handedness] = _new_job(
            workspace, tmp_path / handedness, _shifted_cu3au(handedness), callee, standardize=False
        ).job_id
    _run(workspace)
    results = {handedness: _assert_matches_direct(_succeeded(workspace, job)) for handedness, job in jobs.items()}
    left, right = results["left"]["band_structure"], results["right"]["band_structure"]
    assert left["qpoints"] == right["qpoints"]
    assert numpy.allclose(left["frequencies"], right["frequencies"], atol=1e-9, rtol=0)
    cell = _succeeded(workspace, jobs["left"]).result_file("unitcell/POSCAR").read_text(encoding="utf-8")
    # The left-handed cell is made right-handed by swapping a2 and a3, not by inverting the crystal.
    assert numpy.linalg.det(numpy.array([[float(x) for x in line.split()] for line in cell.splitlines()[2:5]])) > 0
    assert _cartesian_sites(cell) == _cartesian_sites(_shifted_cu3au("left"))


def test_a_failed_displacement_is_recorded_and_fails_the_job(tmp_path: Path) -> None:
    workspace = _workspace(tmp_path)
    marker = tmp_path / "fail-once"
    marker.touch()
    job = _new_job(workspace, tmp_path, CU, _callee(tmp_path / "callee"), call_parameters={"fail_once": str(marker)})
    _run(workspace)
    record = _record(workspace, job.job_id)
    assert record.state == "failed" and record.failure is not None
    assert record.failure.code == "phonons.child_failed" and "displacement-000" in record.failure.message
    [child] = json.loads(record.result_file("phonons.json").read_text(encoding="utf-8"))["children"]
    assert child["label"] == "displacement-000" and child["kind"] == "failed"
    assert child["failure"]["code"] == "tests.forced" and child["job_id"]
