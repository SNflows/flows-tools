"""
This type stub file was generated by pyright.
"""

import abc
from matplotlib.transforms import Transform

__all__ = ['CurvedTransform', 'CoordinateTransform', 'World2PixelTransform', 'Pixel2WorldTransform']
class CurvedTransform(Transform, metaclass=abc.ABCMeta):
    """
    Abstract base class for non-affine curved transforms
    """
    input_dims = ...
    output_dims = ...
    is_separable = ...
    def transform_path(self, path): # -> Path:
        """
        Transform a Matplotlib Path

        Parameters
        ----------
        path : :class:`~matplotlib.path.Path`
            The path to transform

        Returns
        -------
        path : :class:`~matplotlib.path.Path`
            The resulting path
        """
        ...
    
    transform_path_non_affine = ...
    def transform(self, input):
        ...
    
    def inverted(self):
        ...
    


class CoordinateTransform(CurvedTransform):
    has_inverse = ...
    def __init__(self, input_system, output_system) -> None:
        ...
    
    @property
    def same_frames(self):
        ...
    
    @same_frames.setter
    def same_frames(self, same_frames): # -> None:
        ...
    
    def transform(self, input_coords):
        """
        Transform one set of coordinates to another
        """
        ...
    
    transform_non_affine = ...
    def inverted(self): # -> CoordinateTransform:
        """
        Return the inverse of the transform
        """
        ...
    


class World2PixelTransform(CurvedTransform, metaclass=abc.ABCMeta):
    """
    Base transformation from world to pixel coordinates
    """
    has_inverse = ...
    frame_in = ...
    @property
    @abc.abstractmethod
    def input_dims(self): # -> None:
        """
        The number of input world dimensions
        """
        ...
    
    @abc.abstractmethod
    def transform(self, world): # -> None:
        """
        Transform world to pixel coordinates. You should pass in a NxM array
        where N is the number of points to transform, and M is the number of
        dimensions. This then returns the (x, y) pixel coordinates
        as a Nx2 array.
        """
        ...
    
    @abc.abstractmethod
    def inverted(self): # -> None:
        """
        Return the inverse of the transform
        """
        ...
    


class Pixel2WorldTransform(CurvedTransform, metaclass=abc.ABCMeta):
    """
    Base transformation from pixel to world coordinates
    """
    has_inverse = ...
    frame_out = ...
    @property
    @abc.abstractmethod
    def output_dims(self): # -> None:
        """
        The number of output world dimensions
        """
        ...
    
    @abc.abstractmethod
    def transform(self, pixel): # -> None:
        """
        Transform pixel to world coordinates. You should pass in a Nx2 array
        of (x, y) pixel coordinates to transform to world coordinates. This
        will then return an NxM array where M is the number of dimensions.
        """
        ...
    
    @abc.abstractmethod
    def inverted(self): # -> None:
        """
        Return the inverse of the transform
        """
        ...
    


