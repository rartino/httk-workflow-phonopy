# Collecting phonopy calculations

phonopy's finite-displacement method starts from a unit cell, writes displaced
supercells (`phonopy -d --dim 2 2 2` writes `SPOSCAR`, `POSCAR-001`,
`POSCAR-002`, ... and `phonopy_disp.yaml`), and builds force constants from
the forces a calculator computes for each displaced supercell. When those
forces were computed with VASP, one static run per displacement, `httk collect`
recognizes the finished calculation and stores it as one run with the
displacement runs as its children.

## The expected layout

```text
CsCl/
├── phonopy_disp.yaml     # written by phonopy -d (may be compressed)
├── POSCAR                # the unit cell phonopy started from
├── FORCE_SETS            # optional, written by phonopy -f
├── SPOSCAR, POSCAR-001, POSCAR-002, ...   # phonopy's supercells (not read)
├── disp-001/             # a finished VASP static run of POSCAR-001
│   ├── INCAR
│   ├── POSCAR
│   └── OUTCAR
└── disp-002/             # likewise for POSCAR-002
```

The directory is recognized by the `phonopy.calculation.vasp` collector
(priority 20) when it holds `phonopy_disp.yaml` and `POSCAR`, each possibly
compressed, and has exactly one `disp-N` subdirectory per displacement listed
in `phonopy_disp.yaml`, numbered 1 to the number of displacements (`disp-001`
is phonopy's usual naming; the digits may be padded or not). Both phonopy 2
dataset kinds are counted: the usual one (`phonopy -d`, entries
`- atom: N`) and random displacements (`phonopy -d --rd N`, entries `- # N`),
also when a YAML tool has re-indented the list. Otherwise it is
reported as unclaimed with the reason, for example `1 displacement directories
for 2 displacements`. Each `disp-N` directory must be a finished VASP static
run (`NSW = 0` or `IBRION = -1` echoed in its `OUTCAR`), which the
`vasp.calculation.static` collector of *httk-workflow-vasp* recognizes.

Not recognized are phonopy 1.x trees, whose dataset is `disp.yaml`, and trees
made with `phonopy -c POSCAR-unitcell -d`, which have no `POSCAR`; rename the
unit cell to `POSCAR` to collect such a tree.

The calculation's identity is a SHA-256 digest over the digest of
`phonopy_disp.yaml` and `POSCAR` followed, for each displacement directory in
order, by its name and the digest of its `INCAR`, `KPOINTS`, `POTCAR` and
`POSCAR`, the inputs the VASP static collector identifies that run by. All are
digested decompressed, so moving or compressing the directory keeps the
identity, while the same displacements computed with other VASP settings
(another `ENCUT`, k-point mesh or `POTCAR`) are another calculation.

## What is stored

The phonopy directory becomes one run, `phonopy.calculation.vasp:<identity>`,
with:

- the input `initial_structure`: the unit cell read from `POSCAR`, at the VASP
  collectors' precision rule: from the convergence settings an `OUTCAR` beside
  it echoes, and 0.001 Å in the usual case of a phonopy directory without one;
- the output `displacement_dataset`: `phonopy_disp.yaml` as a `files` entry,
  linked as a product of the unit cell;
- the output `force_sets`: `FORCE_SETS` as a `files` entry, when present;
  otherwise the role is reported as unfulfilled, which does not fail the sweep;
- one child-run link per displacement, labelled `disp-001`, `disp-002`, ...,
  to the run of that directory.

A `files` entry records the file's path relative to the swept root as its
`url`, its name, its size and its SHA-256 digest, as stored on disk. The media
type is `application/yaml` for an uncompressed `phonopy_disp.yaml` and unset
otherwise (compressed files and `FORCE_SETS`).

Each displacement run is collected by the VASP static collector as on its own:
its supercell as `initial_structure` and its `total_energy`, stored as a typed
total-energy record. A displacement directory that is missing, not recognized,
or whose own collection degraded (an `OUTCAR` without a final energy, say)
degrades the phonopy run, which then names only the displacement runs that
were collected; the degraded displacement is reported on its own line.

```{admonition} Phonon properties are not stored yet
:class: note

Force constants, phonon band structures, densities of states and
thermodynamic properties have no definitions in the *httk₂* schemas yet, so
they are not computed or stored. What is stored is the provenance of the
calculation, the phonopy dataset files, and the total energy of each
displaced supercell. The forces themselves stay in the `OUTCAR` files and in
`FORCE_SETS`.
```

A phonopy directory that is also a VASP calculation of its own (an `OUTCAR`
beside `POSCAR`, from relaxing the unit cell in place) is claimed by
`phonopy.calculation.vasp`, which has the higher priority; that VASP run is
then not collected, and `--dry-run` names its collector under `also_matched`.
Keep the relaxation in its own directory.

## A sweep

```console
httk collect calculations/ --dry-run
```

prints one `httk-collect-claim` line per recognized directory; the phonopy
directory lists the directories it consumes (identities shortened here):

```text
{"also_matched":[],"collector":"phonopy.calculation.vasp","consumes":["CsCl/disp-001","CsCl/disp-002"],"directory":"CsCl","duplicate_of":null,"format":"httk-collect-claim","format_version":1,"identity":"8f4c5a86…","kind":"claimed","priority":20,"reason":null}
{"also_matched":[],"collector":"vasp.calculation.static","consumes":[],"directory":"CsCl/disp-001","duplicate_of":null,"format":"httk-collect-claim","format_version":1,"identity":"69ff6c78…","kind":"claimed","priority":10,"reason":null}
{"also_matched":[],"collector":"vasp.calculation.static","consumes":[],"directory":"CsCl/disp-002","duplicate_of":null,"format":"httk-collect-claim","format_version":1,"identity":"0746c67d…","kind":"claimed","priority":10,"reason":null}
```

With `disp-002/` missing, the phonopy directory is declined instead:

```text
{"also_matched":[],"collector":"phonopy.calculation.vasp","consumes":[],"directory":"CsCl","duplicate_of":null,"format":"httk-collect-claim","format_version":1,"identity":null,"kind":"unclaimed","priority":null,"reason":"1 displacement directories for 2 displacements"}
```

Collecting into a store collects the displacement runs first and the phonopy
run after them:

```console
httk collect calculations/ --into results.sqlite --id-base mydb
```

The same from Python:

```python
from httk.workflow import collect_tree, store_collected

items = list(collect_tree("calculations"))
reports = store_collected(items, "results.sqlite", id_base="mydb")
```

## Building blocks

`httk.codes.phonopy.collect` has the two helpers the collector's
`recognize.py` hook is made of, for a collector of one's own (a phonopy
calculation computed with another code, say):

- {py:func}`~httk.codes.phonopy.collect.read_displacement_count` counts the
  displacements in a `phonopy_disp.yaml`, reading phonopy's own layout line by
  line, without a YAML parser;
- {py:func}`~httk.codes.phonopy.collect.displacement_directories` lists the
  `disp-N` subdirectories in displacement order.

```python
from pathlib import Path

from httk.codes.phonopy.collect import displacement_directories, read_displacement_count

directory = Path("calculations/CsCl")
print(read_displacement_count(directory / "phonopy_disp.yaml"))  # 2
print(displacement_directories(directory))  # ('disp-001', 'disp-002')
```
