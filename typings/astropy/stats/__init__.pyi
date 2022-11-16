"""
This type stub file was generated by pyright.
"""

from . import bayesian_blocks as _bb, biweight, circstats, funcs, histogram as _hist, info_theory, jackknife, sigma_clipping, spatial
from .funcs import *
from .biweight import *
from .sigma_clipping import *
from .jackknife import *
from .circstats import *
from .bayesian_blocks import *
from .histogram import *
from .info_theory import *
from .spatial import *
from .lombscargle import *
from .bls import *

"""
This subpackage contains statistical tools provided for or used by Astropy.

While the `scipy.stats` package contains a wide range of statistical
tools, it is a general-purpose package, and is missing some that are
particularly useful to astronomy or are used in an atypical way in
astronomy. This package is intended to provide such functionality, but
*not* to replace `scipy.stats` if its implementation satisfies
astronomers' needs.

"""
__all__ = []
