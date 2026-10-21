"""The ``phonopy.calculation.vasp`` collector, run through the walker with the registered collectors."""

import gzip
import hashlib
import shutil
from pathlib import Path
from typing import Any, cast

import pytest

pytest.importorskip("httk.atomistic")

from httk.core import FileRecord
from httk.workflow import DirectoryClaim, claims, collect_tree, store_collected
from httk.workflow.calculations import content_digest

from conftest import ENERGIES, phonopy_disp_yaml, phonopy_tree

PHONOPY = "phonopy.calculation.vasp"
STATIC = "vasp.calculation.static"


def _identity(root: Path) -> str:
    """The documented parent identity: the dataset digest, then each displacement run's VASP inputs."""

    digest = hashlib.sha256(content_digest(root, ["phonopy_disp.yaml", "POSCAR"]).encode() + b"\0")
    for name in ("disp-001", "disp-002"):
        inputs = content_digest(root / name, ["INCAR", "KPOINTS", "POTCAR", "POSCAR"])
        digest.update(name.encode() + b"\0" + inputs.encode() + b"\0")
    return digest.hexdigest()


def _directory(item: Any) -> str:
    return item.record.workdir_path.as_posix()


def test_the_parent_is_claimed_and_consumes_its_displacement_runs(tmp_path: Path) -> None:
    root = phonopy_tree(tmp_path / "CsCl")

    outcomes = {outcome.directory: outcome for outcome in claims(root)}

    assert outcomes["."] == DirectoryClaim(
        ".",
        "claimed",
        PHONOPY,
        20,
        _identity(root),
        consumes=("disp-001", "disp-002"),
    )
    # The parent has no OUTCAR, so the VASP collectors do not match it at all.
    assert outcomes["."].also_matched == ()
    assert [(outcomes[name].kind, outcomes[name].collector) for name in ("disp-001", "disp-002")] == [
        ("claimed", STATIC),
        ("claimed", STATIC),
    ]
    assert set(outcomes) == {".", "disp-001", "disp-002"}


def test_the_children_are_collected_before_their_parent(tmp_path: Path) -> None:
    phonopy_tree(tmp_path / "CsCl")

    first, second, parent = collect_tree(tmp_path)

    assert [(item.record.workspace_id, _directory(item)) for item in (first, second, parent)] == [
        (STATIC, "CsCl/disp-001"),
        (STATIC, "CsCl/disp-002"),
        (PHONOPY, "CsCl"),
    ]
    assert all(item.missing_collector is None for item in (first, second, parent))
    assert [item.outputs["total_energy"].value for item in (first, second)] == pytest.approx(ENERGIES)  # type: ignore[attr-defined]
    assert parent.child_runs == (("disp-001", first.run.source_id), ("disp-002", second.run.source_id))
    assert set(parent.inputs) == {"initial_structure"}
    assert getattr(parent.inputs["initial_structure"], "formula", None) is not None
    assert set(parent.outputs) == {"displacement_dataset"}
    assert parent.unfulfilled == ("force_sets",)
    data = (tmp_path / "CsCl" / "phonopy_disp.yaml").read_bytes()
    assert parent.outputs["displacement_dataset"] == FileRecord(
        "CsCl/phonopy_disp.yaml",
        "phonopy_disp.yaml",
        size=len(data),
        media_type="application/yaml",
        sha256=hashlib.sha256(data).hexdigest(),
    )
    # The dataset is linked as a product of the unit cell it was generated from.
    assert [(link.label, link.source_type, link.target_type) for link in parent.products] == [
        ("displacement_dataset", "structures", "files")
    ]


def test_force_sets_are_kept_when_present(tmp_path: Path) -> None:
    root = phonopy_tree(tmp_path)
    (root / "FORCE_SETS").write_text("16\n2\n\n1\n  0.01 0.0 0.0\n", encoding="utf-8")

    *_children, parent = collect_tree(root)

    assert parent.missing_collector is None and parent.unfulfilled == ()
    force_sets = cast(FileRecord, parent.outputs["force_sets"])
    assert (force_sets.url, force_sets.name, force_sets.media_type) == (
        "FORCE_SETS",
        "FORCE_SETS",
        "application/octet-stream",
    )


def test_the_store_links_the_parent_run_to_the_displacement_runs(tmp_path: Path) -> None:
    store = pytest.importorskip("httk.store")
    from httk.core import Run, TotalEnergyRecord

    phonopy_tree(tmp_path / "tree" / "CsCl")
    items = list(collect_tree(tmp_path / "tree"))
    path = tmp_path / "results.sqlite"

    reports = store_collected(items, str(path), id_base="httk.test")

    assert all("storage_error" not in report for report in reports)
    stored = [cast(dict[str, Any], report["stored"]) for report in reports]
    run_ids = [entry["run"] for entry in stored]
    with store.Backend.sqlite(path) as backend:
        runs = []
        for run_id in run_ids:
            searcher = store.SqlStore(backend).searcher()
            variable = searcher.variable(Run)
            searcher.add(variable.id == run_id)
            runs.append(next(iter(searcher.results(run=variable))).run)
        searcher = store.SqlStore(backend).searcher()
        energies = [row.entry for row in searcher.results(entry=searcher.variable(TotalEnergyRecord))]
        searcher = store.SqlStore(backend).searcher()
        files = [row.entry for row in searcher.results(entry=searcher.variable(FileRecord))]
    *children, parent = runs
    links = sorted((edge.label, edge.entry_id) for edge in parent.artifacts if edge.entry_type == "runs")
    assert links == [("disp-001", run_ids[0]), ("disp-002", run_ids[1])]
    # Each displacement run names its own total energy, stored as a typed record.
    by_id = {energy.id: energy.total_energy for energy in energies}
    named = [{edge.label: edge.entry_id for edge in child.outputs}["total_energy"] for child in children]
    assert [by_id[entry_id] for entry_id in named] == pytest.approx(ENERGIES)
    assert [file.url for file in files] == ["CsCl/phonopy_disp.yaml"]
    assert {edge.label for edge in parent.inputs} == {"initial_structure"}
    assert {edge.label: edge.entry_id for edge in parent.outputs}["displacement_dataset"] == files[0].id


def test_a_missing_displacement_run_leaves_the_parent_unclaimed(tmp_path: Path) -> None:
    root = phonopy_tree(tmp_path)
    shutil.rmtree(root / "disp-002")
    declined: list[DirectoryClaim] = []

    items = list(collect_tree(root, on_unclaimed=declined.append))

    assert [(item.record.workspace_id, _directory(item)) for item in items] == [(STATIC, "disp-001")]
    assert [(outcome.directory, outcome.collector, outcome.reason) for outcome in declined] == [
        (".", PHONOPY, "1 displacement directories for 2 displacements")
    ]


def test_misnumbered_displacement_runs_leave_the_parent_unclaimed(tmp_path: Path) -> None:
    root = phonopy_tree(tmp_path)
    (root / "disp-002").rename(root / "disp-003")

    (parent,) = [outcome for outcome in claims(root) if outcome.directory == "."]

    assert parent.kind == "unclaimed"
    assert parent.reason == "displacement directories disp-001, disp-003 are not numbered 1 to 2"


def test_an_unreadable_dataset_leaves_the_parent_unclaimed(tmp_path: Path) -> None:
    root = phonopy_tree(tmp_path)
    (root / "phonopy_disp.yaml").write_text("phonopy:\n  version: 2\n", encoding="utf-8")

    (parent,) = [outcome for outcome in claims(root) if outcome.directory == "."]

    assert parent.kind == "unclaimed" and "no top-level 'displacements:' list" in (parent.reason or "")


def test_an_unfinished_displacement_run_degrades_it_and_its_parent(tmp_path: Path) -> None:
    phonopy_tree(tmp_path / "CsCl", energies=(ENERGIES[0], None))

    first, second, parent = collect_tree(tmp_path)

    assert first.missing_collector is None
    assert second.missing_collector is not None and second.missing_collector.startswith("CsCl/disp-002: ")
    assert parent.missing_collector is not None
    assert parent.missing_collector.startswith(
        "CsCl: consumed directories not collected: consumed calculation CsCl/disp-002 failed: "
    )
    assert parent.child_runs == (("disp-001", first.run.source_id),)


def test_a_compressed_dataset_is_collected_with_the_same_identity(tmp_path: Path) -> None:
    plain = phonopy_tree(tmp_path / "plain")
    packed = phonopy_tree(tmp_path / "packed")
    (packed / "phonopy_disp.yaml").unlink()
    data = phonopy_disp_yaml().encode()
    with gzip.open(packed / "phonopy_disp.yaml.gz", "wb") as handle:
        handle.write(data)

    (plain_claim,) = [outcome for outcome in claims(plain) if outcome.directory == "."]
    (packed_claim,) = [outcome for outcome in claims(packed) if outcome.directory == "."]
    *_children, parent = collect_tree(packed)

    assert packed_claim.kind == "claimed" and packed_claim.identity == plain_claim.identity
    assert parent.missing_collector is None
    dataset = cast(FileRecord, parent.outputs["displacement_dataset"])
    raw = (packed / "phonopy_disp.yaml.gz").read_bytes()
    assert (dataset.url, dataset.name, dataset.media_type) == (
        "phonopy_disp.yaml.gz",
        "phonopy_disp.yaml.gz",
        "application/octet-stream",
    )
    assert (dataset.size, dataset.sha256) == (len(raw), hashlib.sha256(raw).hexdigest())


def test_the_same_displacements_computed_with_other_vasp_settings_are_another_calculation(tmp_path: Path) -> None:
    phonopy_tree(tmp_path / "encut-520", encut=520)
    phonopy_tree(tmp_path / "encut-600", encut=600)

    parents = {outcome.directory: outcome for outcome in claims(tmp_path) if outcome.collector == PHONOPY}
    items = list(collect_tree(tmp_path))

    assert parents["encut-520"].identity != parents["encut-600"].identity
    assert parents["encut-600"].duplicate_of is None
    collected = [item for item in items if item.record.workspace_id == PHONOPY]
    assert [(_directory(item), item.missing_collector, len(item.child_runs)) for item in collected] == [
        ("encut-520", None, 2),
        ("encut-600", None, 2),
    ]


def test_a_random_displacement_dataset_is_collected(tmp_path: Path) -> None:
    root = phonopy_tree(tmp_path)
    (root / "phonopy_disp.yaml").write_text(phonopy_disp_yaml(random=True), encoding="utf-8")

    *_children, parent = collect_tree(root)

    assert parent.record.workspace_id == PHONOPY and parent.missing_collector is None
    assert [label for label, _source in parent.child_runs] == ["disp-001", "disp-002"]


def test_a_parent_that_is_also_a_vasp_run_is_claimed_by_phonopy(tmp_path: Path) -> None:
    root = phonopy_tree(tmp_path)
    shutil.copy(root / "disp-001" / "OUTCAR", root / "OUTCAR")

    (parent,) = [outcome for outcome in claims(root) if outcome.directory == "."]

    # The documented limitation: the higher priority wins, and the VASP run is not collected.
    assert (parent.kind, parent.collector, parent.also_matched) == ("claimed", PHONOPY, (STATIC,))
