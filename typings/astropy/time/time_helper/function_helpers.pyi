"""
This type stub file was generated by pyright.
"""

import numpy as np

"""
Helpers for overriding numpy functions in
`~astropy.time.Time.__array_function__`.
"""
UNSUPPORTED_FUNCTIONS = ...
CUSTOM_FUNCTIONS = ...
custom_functions = ...
@custom_functions(helps=np.linspace)
def linspace(tstart, tstop, *args, **kwargs): # -> _NotImplementedType | tuple[Any | Unknown, Any | Unknown]:
    ...

