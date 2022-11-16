"""
This type stub file was generated by pyright.
"""

from .core import *
from .quantity import *
from . import astrophys, cgs, misc, photometric, si
from .function import units as function_units
from .si import *
from .astrophys import *
from .photometric import *
from .cgs import *
from .physical import *
from .function.units import *
from .misc import *
from .equivalencies import *
from .function.core import *
from .function.logarithmic import *
from .decorators import *
from .structured import *
def __getattr__(name: str) -> Any: ...
"""
This subpackage contains classes and functions for defining and converting
between different physical units.

This code is adapted from the `pynbody
<https://github.com/pynbody/pynbody>`_ units module written by Andrew
Pontzen, who has granted the Astropy project permission to use the
code under a BSD license.
"""
def __getattr__(attr): # -> Unit | IrreducibleUnit | ((H0: Unknown | None = None) -> Equivalency):
    ...

