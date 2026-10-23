# httk-workflow-phonopy

![Status: Early beta](https://img.shields.io/badge/status-early--beta-orange)

> **⚠️ EARLY BETA**
>
> This is an early beta release of *httk₂*. The organization of the packages
> and their APIs should not yet be regarded as stable, and may change between
> releases.

*httk-workflow-phonopy* adds [phonopy](https://phonopy.github.io/phonopy/)
support to [*httk-workflow*](https://github.com/httk/httk-workflow), the
workflow engine of [*httk₂*](https://github.com/httk/httk2), in two parts:

- the workflow package `phonopy.phonons` (in `workflows/phonopy-phonons`),
  which computes harmonic phonons by finite displacements with forces from any
  force workflow the job binds (see [phonopy.phonons](#phonopyphonons));
- `httk.codes.phonopy`, which registers the `phonopy.calculation.vasp`
  collector so that `httk collect` recognizes finished phonopy
  finite-displacement calculations whose displaced supercells were computed
  with VASP; nothing needs to be configured.

phonopy itself is not included: `phonopy.phonons` declares it in its
`requires`, and the collector only reads phonopy's files.

A phonopy calculation directory holds `phonopy_disp.yaml`, the unit cell
`POSCAR`, optionally `FORCE_SETS`, and one finished VASP static run per
displacement in `disp-001/`, `disp-002/`, and so on. It is collected as one run
with the unit cell as its input and the phonopy dataset files as `files`
outputs, linked to the displacement runs as its child runs; each displacement
run is collected by the VASP static collector of
[*httk-workflow-vasp*](https://github.com/httk/httk-workflow-vasp), with its
total energy. The collector stores no phonon properties (force constants,
band structures); `phonopy.phonons` computes and records them for its own runs.

## Install

```console
python -m pip install httk-workflow-phonopy
```

Reading structures needs *httk-atomistic* as well, as for the VASP collectors.

## Collecting phonopy calculations

```console
httk collect calculations/ --dry-run
httk collect calculations/ --into results.sqlite --id-base mydb
```

The expected layout and what is stored are documented in
[`docs/usage.md`](docs/usage.md) and at
[docs.httk.org/httk-workflow-phonopy](https://docs.httk.org/httk-workflow-phonopy/).

## phonopy.phonons

`phonopy.phonons` computes the harmonic phonons of one structure by finite
displacements with phonopy. It computes no forces itself: its call slot
`forces` takes any workflow that turns a structure into the forces on its sites
(the interface `workflows/phonopy-phonons/calls/structure-forces.json`: input
`structure`, output `forces`, a `cartesian_site_forces` record that is a
product of that structure). By default it is `vasp.static` of
*workflows-vasp*, configured with the settings `vasp.phonons` gives its
displacement runs; a job can bind `ase.static` of *workflows-examples* or any
other fitting workflow instead.

Scope:

- The input `structure` must already be relaxed, with the same force model and
  settings as the binding. There is no relaxation and no non-analytical
  correction (LO-TO splitting is not modelled for polar materials).
- `displace` standardizes the structure to its HPKOT primitive cell (as
  `httk.atomistic.band_path` gives it; `standardize = false` keeps the given
  cell), makes it right-handed by swapping `a2` and `a3` (the same crystal, not
  its mirror image; `vasp.phonons` still inverts the crystal instead), and
  orders its sites by species with a stable sort. That cell is the
  `phonon_cell` output. phonopy builds a diagonal supercell of perpendicular
  widths of at least `min_length` and its displacements. Each displaced
  supercell is written, read back and checked against phonopy's site order
  (`phonons.input_invalid` otherwise), then called through the slot as child
  `displacement-NNN`.
- `analyse` waits for every displacement to finish. A failed one
  (`phonons.child_failed`) or one without a child job
  (`phonons.child_missing`) is recorded in `phonons.json` with its job id and
  failure before the job fails. It then reads each child's forces, refuses
  unreadable outputs (`phonons.forces_unavailable`) or forces that are not one
  finite 3-vector per supercell atom (`phonons.forces_invalid`), and computes
  the force constants, the band structure along the HPKOT path, the density of
  states and the harmonic thermodynamics. Every force array used is kept in
  `phonons.json`, so a retried analysis reads the same values.

Install the plugin from this repository (`httk plugin install .`) or use the
package directly with `--workflow-dir workflows/phonopy-phonons`.

### Binding the forces

By default the slot is bound to `vasp.static`, which must be installed from
*workflows-vasp* (`httk plugin install` of that repository, or
`httk workflow install vasp.static`, or `job new --install`). Installing
`phonopy.phonons` does not install it; it only prints a notice when it is
absent. The binding passes these `vasp.static` parameters, the ones
`vasp.phonons` applies to its displacement runs:

| Parameter | Value | Why |
| --- | --- | --- |
| `excluded_remedies` | `bump_kpoints equal_kpoints centering incar.ISMEAR incar.SIGMA incar.LREAL scale_lattice` | These error remedies would change the k-point mesh, its centering, the smearing, the projection scheme or the cell of one child but not the others, and force constants built from forces of different numerics are meaningless. A child that would need one fails instead. |
| `kpoint_density` | `20.0` | One fixed density, so every supercell (all have the same lattice) gets the same automatic mesh. |
| `centering` | `Gamma` | The mesh centering `vasp.phonons` uses. |
| `incar_tags` | `{"EDIFF": "1E-7", "LREAL": ".FALSE.", "ADDGRID": ".TRUE.", "ISMEAR": 0, "SIGMA": 0.05, "ICHARG": 2}` | Where `vasp.phonons`' `INCAR.phonon` differs from `vasp.static`'s `INCAR.static`: tight `EDIFF`, no real-space projection, the additional grid, Gaussian smearing, a start from atomic charges. |

```console
httk job new --install --workflow phonopy.phonons --input structure=POSCAR
httk workflow run
httk collect --into results.sqlite
```

`--call-parameter forces NAME=VALUE` overrides one entry and keeps the others,
for example `--call-parameter forces kpoint_density=30`; other `vasp.static`
parameters (`accuracy_per_atom`, `pseudopotential_library`, ...) are added the
same way. An overridden `incar_tags` replaces the whole table, so repeat the
tags above that you want to keep.

Bind another workflow with `--call forces=WORKFLOW`; the parameters above are
then not passed. With ASE's EMT potential (no external code; useful to try the
workflow), from the *workflows-examples* repository (`httk plugin install` of
that repository or `httk workflow install`):

```console
httk job new --install --workflow phonopy.phonons --input structure=POSCAR \
    --call forces=ase.static --call-parameter forces calculator=ase.calculators.emt:EMT
```

Code-specific numerical consistency across the displacement children is the
binding's responsibility: a code-neutral workflow cannot judge one code's
settings, so `phonopy.phonons` does not compare the children's inputs (unlike
`vasp.phonons`) and records no remedy histories, whose location is
code-specific. Each child stays a job of its own workflow, with its own inputs
and records.

### Parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `standardize` | `true` | Use the HPKOT primitive cell; `false` keeps the given cell, which must be a basis of its HPKOT lattice. |
| `min_length` | `20.0` | Minimum perpendicular width of the diagonal supercell, in Å. |
| `max_atoms` | `300` | Refuse (`phonons.supercell_too_large`) a larger supercell. |
| `displacement` | `0.01` | phonopy's displacement amplitude, in Å. |
| `symmetry_tolerance` | `1e-3` | Symmetry tolerance of the HPKOT cell and path, in Å. |
| `mesh_length` | `100.0` | phonopy's q-point mesh length for the DOS and the thermodynamics. |
| `band_points` | `51` | q-points per segment of the band path. |
| `temperatures` | 0 to 1000 K in steps of 10 K | Temperatures of the thermodynamics. |
| `cutoff_frequency` | `0.05` | Modes below this frequency (THz) are left out of the thermodynamics. |
| `imaginary` | `raise` | `raise`: no thermodynamics with imaginary modes; `omit`: leave them out. |
| `imaginary_tolerance` | `0.05` | Mesh frequencies below minus this (THz) count as imaginary. |

The names and defaults are those of `vasp.phonons`.

### Outputs

- `phonon_cell` (structure): the species-ordered HPKOT cell the phonon records
  refer to;
- `phonon_band_structure`, `phonon_density_of_states`,
  `vibrational_thermodynamics` and `zero_point_energy` (records of the phonon
  cell; the thermodynamic ones only when computed);
- `phonons_result` (`phonons.json`: the forces of every child, the supercell,
  the settings and the results), `phonopy_disp` (`phonopy_disp.yaml`) and
  `force_sets` (`FORCE_SETS`) as files.

Collecting needs *httk-atomistic* and *httk-analyse*; the displacement
children are collected as child runs by their own workflow.

## Running tests

`make test` runs the normal profile; `make ci` runs formatting, lint, both type
checkers and the extended tests. The collector tests build a synthetic phonopy
tree with VASP static runs and need neither phonopy nor VASP. The
`phonopy.phonons` tests run real jobs with an ASE EMT callee and skip without
phonopy and ASE; when the *workflows-examples* checkout sits beside this
repository, they also run against its `ase.static`.
