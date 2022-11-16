"""
This type stub file was generated by pyright.
"""

from .nddata import NDData
from .mixins.ndslicing import NDSlicingMixin
from .mixins.ndarithmetic import NDArithmeticMixin
from .mixins.ndio import NDIOMixin

__all__ = ['NDDataArray']
class NDDataArray(NDArithmeticMixin, NDSlicingMixin, NDIOMixin, NDData):
    """
    An ``NDData`` object with arithmetic. This class is functionally equivalent
    to ``NDData`` in astropy  versions prior to 1.0.

    The key distinction from raw numpy arrays is the presence of
    additional metadata such as uncertainties, a mask, units, flags,
    and/or a coordinate system.

    See also: https://docs.astropy.org/en/stable/nddata/

    Parameters
    ----------
    data : ndarray or `NDData`
        The actual data contained in this `NDData` object. Not that this
        will always be copies by *reference* , so you should make copy
        the ``data`` before passing it in if that's the  desired behavior.

    uncertainty : `~astropy.nddata.NDUncertainty`, optional
        Uncertainties on the data.

    mask : array-like, optional
        Mask for the data, given as a boolean Numpy array or any object that
        can be converted to a boolean Numpy array with a shape
        matching that of the data. The values must be ``False`` where
        the data is *valid* and ``True`` when it is not (like Numpy
        masked arrays). If ``data`` is a numpy masked array, providing
        ``mask`` here will causes the mask from the masked array to be
        ignored.

    flags : array-like or `~astropy.nddata.FlagCollection`, optional
        Flags giving information about each pixel. These can be specified
        either as a Numpy array of any type (or an object which can be converted
        to a Numpy array) with a shape matching that of the
        data, or as a `~astropy.nddata.FlagCollection` instance which has a
        shape matching that of the data.

    wcs : None, optional
        WCS-object containing the world coordinate system for the data.

        .. warning::
            This is not yet defined because the discussion of how best to
            represent this class's WCS system generically is still under
            consideration. For now just leave it as None

    meta : `dict`-like object, optional
        Metadata for this object.  "Metadata" here means all information that
        is included with this object but not part of any other attribute
        of this particular object.  e.g., creation date, unique identifier,
        simulation parameters, exposure time, telescope name, etc.

    unit : `~astropy.units.UnitBase` instance or str, optional
        The units of the data.


    Raises
    ------
    ValueError :
        If the `uncertainty` or `mask` inputs cannot be broadcast (e.g., match
        shape) onto ``data``.
    """
    def __init__(self, data, *args, flags=..., **kwargs) -> None:
        ...
    
    @property
    def uncertainty(self): # -> NDUncertainty:
        ...
    
    @uncertainty.setter
    def uncertainty(self, value): # -> None:
        ...
    
    @property
    def unit(self): # -> Unit | None:
        ...
    
    @unit.setter
    def unit(self, value): # -> None:
        ...
    
    @property
    def mask(self): # -> NDArray[bool_] | bool_ | None:
        ...
    
    @mask.setter
    def mask(self, value): # -> None:
        ...
    
    @property
    def shape(self): # -> Any:
        """
        shape tuple of this object's data.
        """
        ...
    
    @property
    def size(self): # -> Any:
        """
        integer size of this object's data.
        """
        ...
    
    @property
    def dtype(self): # -> Any:
        """
        `numpy.dtype` of this object's data.
        """
        ...
    
    @property
    def ndim(self): # -> Any:
        """
        integer dimensions of this object's data
        """
        ...
    
    @property
    def flags(self): # -> FlagCollection | NDArray[Unknown]:
        ...
    
    @flags.setter
    def flags(self, value): # -> None:
        ...
    
    def __array__(self): # -> masked_array | NDArray[Any]:
        """
        This allows code that requests a Numpy array to use an NDData
        object as a Numpy array.
        """
        ...
    
    def __array_prepare__(self, array, context=...): # -> masked_array:
        """
        This ensures that a masked array is returned if self is masked.
        """
        ...
    
    def convert_unit_to(self, unit, equivalencies=...): # -> Self@NDDataArray:
        """
        Returns a new `NDData` object whose values have been converted
        to a new unit.

        Parameters
        ----------
        unit : `astropy.units.UnitBase` instance or str
            The unit to convert to.

        equivalencies : list of tuple
           A list of equivalence pairs to try if the units are not
           directly convertible.  See :ref:`astropy:unit_equivalencies`.

        Returns
        -------
        result : `~astropy.nddata.NDData`
            The resulting dataset

        Raises
        ------
        `~astropy.units.UnitsError`
            If units are inconsistent.

        """
        ...
    


