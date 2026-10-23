# *httk-workflow-phonopy*

This site documents the *httk-workflow-phonopy* module. For the full
documentation of *httk₂*, see [docs.httk.org](https://docs.httk.org).

The module ships the workflow package `phonopy.phonons`, code-neutral
harmonic phonons by finite displacements whose forces come from any force
workflow the job binds (see {doc}`usage`), and teaches `httk collect` to
recognize finished
[phonopy](https://phonopy.github.io/phonopy/) finite-displacement
calculations whose displaced supercells were computed with VASP. Such a
calculation is a directory holding phonopy's `phonopy_disp.yaml` and the unit
cell `POSCAR`, with one VASP static run per displacement in `disp-001/`,
`disp-002/`, and so on. It is collected as one run that links the displacement
runs as its child runs; the displacement runs themselves are collected by the
`vasp.calculation.static` collector of *httk-workflow-vasp*.

Installing it registers the `phonopy.calculation.vasp` collector through the
`httk.registry.codes.phonopy` registration package. It registers no code: the
collector only reads phonopy's files. The building blocks of the
collector are in `httk.codes.phonopy.collect`. The collector stores no phonon
properties (force constants, band structures, densities of states);
`phonopy.phonons` records them for its own runs; see {doc}`usage`.

```{admonition} Quick links
:class: tip

- {doc}`usage` — the expected layout, what is stored, a sweep, and the `phonopy.phonons` workflow
- {doc}`reference/index` — the generated API reference
```

## Install

```console
python -m pip install httk-workflow-phonopy
```

Reading the structures needs *httk-atomistic* as well, as for the VASP
collectors.

```{toctree}
:maxdepth: 2
:caption: Documentation

usage
reference/index
```
