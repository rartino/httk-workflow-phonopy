"""Collect hook for ``phonopy.calculation.vasp``.

The unit cell is read from the directory's POSCAR, and the phonopy dataset
files are kept as ``files`` entries. The displacement runs are collected by the
VASP static collector and linked as child runs by the framework.
"""

from pathlib import Path, PurePosixPath

from httk.codes.vasp.collect import read_structure, structure_precision
from httk.core import FileRecord
from httk.core.datastream.compression import split_compression_suffix
from httk.core.digests import sha256_file
from httk.workflow.collecting import JobRecord, existing_file


def _file(where: PurePosixPath, path: Path, media_type: str | None) -> FileRecord:
    """Describe one dataset file of the directory *where*, relative to the swept root."""

    compressed = split_compression_suffix(path.name)[1] is not None
    return FileRecord(
        (where / path.name).as_posix(),
        path.name,
        size=path.stat().st_size,
        media_type="application/octet-stream" if compressed or media_type is None else media_type,
        sha256=sha256_file(path),
    )


def collect(record: JobRecord) -> dict[str, object]:
    """Return the unit cell and the displacement dataset, with ``FORCE_SETS`` when present.

    :param record: The stand-in job record of the directory.
    :return: The ``initial_structure`` input role and the ``displacement_dataset``
        and, when present, ``force_sets`` output roles.
    """

    where = record.workdir_path
    if where is None:
        raise ValueError("the job has no working directory")
    poscar = record.result_file("POSCAR")
    roles: dict[str, object] = {
        "initial_structure": read_structure(poscar, precision=structure_precision(poscar.parent)),
        "displacement_dataset": _file(where, record.result_file("phonopy_disp.yaml"), "application/yaml"),
    }
    force_sets = existing_file(poscar.parent / "FORCE_SETS")
    if force_sets is not None:
        roles["force_sets"] = _file(where, force_sets, None)
    return roles
