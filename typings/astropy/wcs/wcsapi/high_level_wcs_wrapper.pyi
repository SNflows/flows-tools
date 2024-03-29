"""
This type stub file was generated by pyright.
"""

from .high_level_api import HighLevelWCSMixin

__all__ = ['HighLevelWCSWrapper']
class HighLevelWCSWrapper(HighLevelWCSMixin):
    """
    Wrapper class that can take any :class:`~astropy.wcs.wcsapi.BaseLowLevelWCS`
    object and expose the high-level WCS API.
    """
    def __init__(self, low_level_wcs) -> None:
        ...
    
    @property
    def low_level_wcs(self): # -> BaseLowLevelWCS:
        ...
    
    @property
    def pixel_n_dim(self): # -> None:
        """
        See `~astropy.wcs.wcsapi.BaseLowLevelWCS.world_n_dim`
        """
        ...
    
    @property
    def world_n_dim(self): # -> None:
        """
        See `~astropy.wcs.wcsapi.BaseLowLevelWCS.world_n_dim`
        """
        ...
    
    @property
    def world_axis_physical_types(self): # -> None:
        """
        See `~astropy.wcs.wcsapi.BaseLowLevelWCS.world_axis_physical_types`
        """
        ...
    
    @property
    def world_axis_units(self): # -> None:
        """
        See `~astropy.wcs.wcsapi.BaseLowLevelWCS.world_axis_units`
        """
        ...
    
    @property
    def array_shape(self): # -> None:
        """
        See `~astropy.wcs.wcsapi.BaseLowLevelWCS.array_shape`
        """
        ...
    
    @property
    def pixel_bounds(self): # -> None:
        """
        See `~astropy.wcs.wcsapi.BaseLowLevelWCS.pixel_bounds`
        """
        ...
    
    @property
    def axis_correlation_matrix(self):
        """
        See `~astropy.wcs.wcsapi.BaseLowLevelWCS.axis_correlation_matrix`
        """
        ...
    
    def __str__(self) -> str:
        ...
    
    def __repr__(self): # -> str:
        ...
    


