"""Register the phonopy calculation collector implemented by :mod:`httk.codes.phonopy`.

phonopy is only collected, never run, so no code (bridge or Bash API) is registered.
"""

from httk.core.register import register_collector

register_collector("phonopy.calculation.vasp", package="httk.codes.phonopy:collectors/phonopy-vasp")
