# httk-workflow-phonopy

![Status: Early beta](https://img.shields.io/badge/status-early--beta-orange)

> **⚠️ EARLY BETA**
>
> This is an early beta release of *httk₂*. The organization of the packages
> and their APIs should not yet be regarded as stable, and may change between
> releases.

*httk-workflow-phonopy* teaches `httk collect` of
[*httk-workflow*](https://github.com/httk/httk-workflow), the workflow engine
of [*httk₂*](https://github.com/httk/httk2), to recognize finished
[phonopy](https://phonopy.github.io/phonopy/) finite-displacement calculations
whose displaced supercells were computed with VASP. It provides
`httk.codes.phonopy` and registers the `phonopy.calculation.vasp` collector;
nothing needs to be configured. phonopy itself is not included, and nothing is
run: the package only collects.

A phonopy calculation directory holds `phonopy_disp.yaml`, the unit cell
`POSCAR`, optionally `FORCE_SETS`, and one finished VASP static run per
displacement in `disp-001/`, `disp-002/`, and so on. It is collected as one run
with the unit cell as its input and the phonopy dataset files as `files`
outputs, linked to the displacement runs as its child runs; each displacement
run is collected by the VASP static collector of
[*httk-workflow-vasp*](https://github.com/httk/httk-workflow-vasp), with its
total energy. Phonon properties (force constants, band structures) are not
stored yet, because the *httk₂* schemas have no definitions for them.

## Install

```console
python -m pip install httk-workflow-phonopy
```

Reading structures needs *httk-atomistic* as well, as for the VASP collectors.

## Use

```console
httk collect calculations/ --dry-run
httk collect calculations/ --into results.sqlite --id-base mydb
```

The expected layout and what is stored are documented in
[`docs/usage.md`](docs/usage.md) and at
[docs.httk.org/httk-workflow-phonopy](https://docs.httk.org/httk-workflow-phonopy/).

## Running tests

`make test` runs the normal profile; `make ci` runs formatting, lint, both type
checkers and the extended tests. The tests build a synthetic phonopy tree with
VASP static runs; they need neither phonopy nor VASP.
