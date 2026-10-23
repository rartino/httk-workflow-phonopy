"""Collect hook for the packaged ``phonopy.phonons`` workflow.

``analyse`` writes ``phonons.json`` (every input of the analysis and its
results), ``phonopy_disp.yaml`` and ``FORCE_SETS``, and ``displace`` the phonon
unit cell ``unitcell/POSCAR``, into the workdir. The records are rebuilt from
the stored values, never recomputed, so collecting twice gives the same content
ids. The displacement children collect as child runs of their own, by their own
workflow's hook.
"""

import json
from pathlib import Path

import httk.core
from httk.core import FileRecord
from httk.core.digests import sha256_file

_MEDIA_TYPES = {"phonons.json": "application/json", "phonopy_disp.yaml": "application/yaml"}


def _file(record, path: Path) -> FileRecord:
    """Describe one result file, by its path relative to the workspace root."""

    return FileRecord(
        path.resolve().relative_to(record.workspace_root.resolve()).as_posix(),
        path.name,
        size=path.stat().st_size,
        media_type=_MEDIA_TYPES.get(path.name, "application/octet-stream"),
        sha256=sha256_file(path),
    )


def collect(record):
    """Return the phonon cell, its band structure, DOS and, when computed, thermodynamics, and the result files.

    :param record: The collected job record.
    :return: Extracted output roles.
    :raises ValueError: If httk-atomistic or httk-analyse is not installed, or ``phonons.json`` holds no phonon
        results.
    """
    try:
        from httk.analyse import definitions, property_records
        from httk.analyse.matsci.phonons import HarmonicThermodynamics
        from httk.analyse.records import bound_values, records
        from httk.atomistic import UnitcellStructureView
    except ImportError as exc:
        # Without httk-analyse the values would be stored as untyped rows a store must later rebuild.
        raise ValueError(
            "collecting phonopy.phonons needs httk-atomistic and httk-analyse; pip install them where you collect"
        ) from exc
    path = record.result_file("phonons.json")
    document = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(document.get("band_structure"), dict) or not isinstance(document.get("dos"), dict):
        raise ValueError(f"{path}: no phonon results")
    kinds = property_records.RECORD_KINDS
    cell = record.result_file("unitcell/POSCAR")
    outputs = {
        # phonopy writes 16 decimals; the cell is stored as written.
        "phonon_cell": UnitcellStructureView(httk.core.load(str(cell), format="vasp-poscar", precision=1e-9)),
        "phonon_band_structure": kinds[definitions.PHONON_BAND_STRUCTURE].from_value(document["band_structure"]),
        "phonon_density_of_states": kinds[definitions.PHONON_DENSITY_OF_STATES].from_value(document["dos"]),
        "phonons_result": _file(record, path),
        "phonopy_disp": _file(record, record.result_file("phonopy_disp.yaml")),
        "force_sets": _file(record, record.result_file("FORCE_SETS")),
    }
    stored = document.get("thermodynamics")
    # Omitted (imaginary modes, or none left): the roles stay unfulfilled, never null records.
    if stored is not None:
        thermodynamics = HarmonicThermodynamics(**stored)
        for bound, value in zip(bound_values(thermodynamics), records(thermodynamics), strict=True):
            outputs[bound.binding.definition.rsplit("/", 1)[1]] = value
    return outputs
