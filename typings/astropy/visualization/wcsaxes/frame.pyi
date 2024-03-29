"""
This type stub file was generated by pyright.
"""

import abc
from collections import OrderedDict

__all__ = ['RectangularFrame1D', 'Spine', 'BaseFrame', 'RectangularFrame', 'EllipticalFrame']
class Spine:
    """
    A single side of an axes.

    This does not need to be a straight line, but represents a 'side' when
    determining which part of the frame to put labels and ticks on.
    """
    def __init__(self, parent_axes, transform) -> None:
        ...
    
    @property
    def data(self): # -> None:
        ...
    
    @data.setter
    def data(self, value): # -> None:
        ...
    
    @property
    def pixel(self): # -> None:
        ...
    
    @pixel.setter
    def pixel(self, value): # -> None:
        ...
    
    @property
    def world(self): # -> None:
        ...
    
    @world.setter
    def world(self, value): # -> None:
        ...
    


class SpineXAligned(Spine):
    """
    A single side of an axes, aligned with the X data axis.

    This does not need to be a straight line, but represents a 'side' when
    determining which part of the frame to put labels and ticks on.
    """
    @property
    def data(self): # -> None:
        ...
    
    @data.setter
    def data(self, value): # -> None:
        ...
    
    @property
    def pixel(self): # -> None:
        ...
    
    @pixel.setter
    def pixel(self, value): # -> None:
        ...
    


class BaseFrame(OrderedDict, metaclass=abc.ABCMeta):
    """
    Base class for frames, which are collections of
    :class:`~astropy.visualization.wcsaxes.frame.Spine` instances.
    """
    spine_class = Spine
    def __init__(self, parent_axes, transform, path=...) -> None:
        ...
    
    @property
    def origin(self): # -> Literal['lower', 'upper']:
        ...
    
    @property
    def transform(self): # -> Unknown:
        ...
    
    @transform.setter
    def transform(self, value): # -> None:
        ...
    
    @property
    def patch(self): # -> PathPatch:
        ...
    
    def draw(self, renderer): # -> None:
        ...
    
    def sample(self, n_samples): # -> OrderedDict[Unknown, Unknown]:
        ...
    
    def set_color(self, color): # -> None:
        """
        Sets the color of the frame.

        Parameters
        ----------
        color : str
            The color of the frame.
        """
        ...
    
    def get_color(self):
        ...
    
    def set_linewidth(self, linewidth): # -> None:
        """
        Sets the linewidth of the frame.

        Parameters
        ----------
        linewidth : float
            The linewidth of the frame in points.
        """
        ...
    
    def get_linewidth(self):
        ...
    
    @abc.abstractmethod
    def update_spines(self):
        ...
    


class RectangularFrame1D(BaseFrame):
    """
    A classic rectangular frame.
    """
    spine_names = ...
    spine_class = SpineXAligned
    def update_spines(self): # -> None:
        ...
    
    def draw(self, renderer): # -> None:
        ...
    


class RectangularFrame(BaseFrame):
    """
    A classic rectangular frame.
    """
    spine_names = ...
    def update_spines(self): # -> None:
        ...
    


class EllipticalFrame(BaseFrame):
    """
    An elliptical frame.
    """
    spine_names = ...
    def update_spines(self): # -> None:
        ...
    
    def draw(self, renderer): # -> None:
        """Override to draw only the outer ellipse,
        not the major and minor axes in the middle.

        FIXME: we may want to add a general method to give the user control
        over which spines are drawn."""
        ...
    


