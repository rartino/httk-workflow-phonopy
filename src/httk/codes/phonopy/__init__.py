"""phonopy calculation support for *httk₂* workflows: the *httk-workflow-phonopy* package.

phonopy computes phonons from the forces of displaced supercells. This package
teaches ``httk collect`` to recognize a finished finite-displacement phonopy
calculation whose displaced supercells were computed with VASP: the
``phonopy.calculation.vasp`` collector, shipped as package data and registered
through ``httk.registry.codes.phonopy``. ``collect`` holds the building blocks
it is made of. phonopy has no runner support here: nothing is run, only
collected. phonopy itself is not included.
"""

from httk.core import register_citation

register_citation(
    applies_to="Phonon calculations with phonopy",
    references=(
        {
            "authors": (
                {"name": "Atsushi Togo"},
                {"name": "Laurent Chaput"},
                {"name": "Terumasa Tadano"},
                {"name": "Isao Tanaka"},
            ),
            "title": "Implementation strategies in phonopy and phono3py",
            "journal": "Journal of Physics: Condensed Matter",
            "volume": "35",
            "pages": "353001",
            "year": "2023",
            "doi": "10.1088/1361-648X/acd831",
            "bib_type": "article",
        },
        {
            "authors": ({"name": "Atsushi Togo"}, {"name": "Isao Tanaka"}),
            "title": "First principles phonon calculations in materials science",
            "journal": "Scripta Materialia",
            "volume": "108",
            "pages": "1-5",
            "year": "2015",
            "doi": "10.1016/j.scriptamat.2015.07.021",
            "bib_type": "article",
        },
    ),
)

__all__: list[str] = []
