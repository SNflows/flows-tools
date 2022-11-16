"""
This type stub file was generated by pyright.
"""

from astropy.utils.decorators import lazyproperty
from .base import BaseWCSWrapper

__all__ = ['sanitize_slices', 'SlicedLowLevelWCS']
def sanitize_slices(slices, ndim): # -> list[Unknown]:
    """
    Given a slice as input sanitise it to an easier to parse format.format

    This function returns a list ``ndim`` long containing slice objects (or ints).
    """
    ...

def combine_slices(slice1, slice2): # -> Integral | Any | slice:
    """
    Given two slices that can be applied to a 1-d array, find the resulting
    slice that corresponds to the combination of both slices. We assume that
    slice2 can be an integer, but slice1 cannot.
    """
    ...

class SlicedLowLevelWCS(BaseWCSWrapper):
    """
    A Low Level WCS wrapper which applies an array slice to a WCS.

    This class does not modify the underlying WCS object and can therefore drop
    coupled dimensions as it stores which pixel and world dimensions have been
    sliced out (or modified) in the underlying WCS and returns the modified
    results on all the Low Level WCS methods.

    Parameters
    ----------
    wcs : `~astropy.wcs.wcsapi.BaseLowLevelWCS`
        The WCS to slice.
    slices : `slice` or `tuple` or `int`
        A valid array slice to apply to the WCS.

    """
    def __init__(self, wcs, slices) -> None:
        ...
    
    @lazyproperty
    def dropped_world_dimensions(self): # -> dict[Unknown, list[Unknown]]:
        """
        Information describing the dropped world dimensions.
        """
        ...
    
    @property
    def pixel_n_dim(self): # -> int:
        ...
    
    @property
    def world_n_dim(self): # -> int:
        ...
    
    @property
    def world_axis_physical_types(self): # -> list[Unknown]:
        ...
    
    @property
    def world_axis_units(self): # -> list[Unknown]:
        ...
    
    @property
    def pixel_axis_names(self): # -> list[Unknown]:
        ...
    
    @property
    def world_axis_names(self): # -> list[Unknown]:
        ...
    
    def pixel_to_world_values(self, *pixel_arrays): # -> ndarray[Unknown, Unknown] | Any | list[Unknown | Any]:
        ...
    
    def world_to_pixel_values(self, *world_arrays): # -> <subclass of list and ndarray> | tuple[Unknown, ...]:
        ...
    
    @property
    def world_axis_object_components(self): # -> list[Unknown]:
        ...
    
    @property
    def world_axis_object_classes(self): # -> dict[Unknown, Unknown]:
        ...
    
    @property
    def array_shape(self): # -> Any | None:
        ...
    
    @property
    def pixel_shape(self): # -> tuple[Any, ...] | None:
        ...
    
    @property
    def pixel_bounds(self): # -> tuple[Unknown, ...] | None:
        ...
    
    @property
    def axis_correlation_matrix(self):
        ...
    

