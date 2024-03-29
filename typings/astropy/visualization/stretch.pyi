"""
This type stub file was generated by pyright.
"""

from .transform import BaseTransform, CompositeTransform

"""
Classes that deal with stretching, i.e. mapping a range of [0:1] values onto
another set of [0:1] values with a transformation
"""
__all__ = ["BaseStretch", "LinearStretch", "SqrtStretch", "PowerStretch", "PowerDistStretch", "SquaredStretch", "LogStretch", "AsinhStretch", "SinhStretch", "HistEqStretch", "ContrastBiasStretch", "CompositeStretch"]
class BaseStretch(BaseTransform):
    """
    Base class for the stretch classes, which, when called with an array
    of values in the range [0:1], return an transformed array of values,
    also in the range [0:1].
    """
    def __add__(self, other): # -> CompositeStretch:
        ...
    
    def __call__(self, values, clip=..., out=...): # -> None:
        """
        Transform values using this stretch.

        Parameters
        ----------
        values : array-like
            The input values, which should already be normalized to the
            [0:1] range.
        clip : bool, optional
            If `True` (default), values outside the [0:1] range are
            clipped to the [0:1] range.
        out : ndarray, optional
            If specified, the output values will be placed in this array
            (typically used for in-place calculations).

        Returns
        -------
        result : ndarray
            The transformed values.
        """
        ...
    
    @property
    def inverse(self): # -> None:
        """A stretch object that performs the inverse operation."""
        ...
    


class LinearStretch(BaseStretch):
    """
    A linear stretch with a slope and offset.

    The stretch is given by:

    .. math::
        y = slope x + intercept

    Parameters
    ----------
    slope : float, optional
        The ``slope`` parameter used in the above formula.  Default is 1.
    intercept : float, optional
        The ``intercept`` parameter used in the above formula.  Default is 0.
    """
    def __init__(self, slope=..., intercept=...) -> None:
        ...
    
    def __call__(self, values, clip=..., out=...): # -> NDArray[Unknown]:
        ...
    
    @property
    def inverse(self): # -> LinearStretch:
        """A stretch object that performs the inverse operation."""
        ...
    


class SqrtStretch(BaseStretch):
    r"""
    A square root stretch.

    The stretch is given by:

    .. math::
        y = \sqrt{x}
    """
    def __call__(self, values, clip=..., out=..., invalid=...): # -> NDArray[Unknown]:
        """
        Transform values using this stretch.

        Parameters
        ----------
        values : array-like
            The input values, which should already be normalized to the
            [0:1] range.
        clip : bool, optional
            If `True` (default), values outside the [0:1] range are
            clipped to the [0:1] range.
        out : ndarray, optional
            If specified, the output values will be placed in this array
            (typically used for in-place calculations).
        invalid : None or float, optional
            Value to assign NaN values generated by this class.  NaNs in
            the input ``values`` array are not changed.  This option is
            generally used with matplotlib normalization classes, where
            the ``invalid`` value should map to the matplotlib colormap
            "under" value (i.e., any finite value < 0).  If `None`, then
            NaN values are not replaced.  This keyword has no effect if
            ``clip=True``.

        Returns
        -------
        result : ndarray
            The transformed values.
        """
        ...
    
    @property
    def inverse(self): # -> PowerStretch:
        """A stretch object that performs the inverse operation."""
        ...
    


class PowerStretch(BaseStretch):
    r"""
    A power stretch.

    The stretch is given by:

    .. math::
        y = x^a

    Parameters
    ----------
    a : float
        The power index (see the above formula).  ``a`` must be greater
        than 0.
    """
    def __init__(self, a) -> None:
        ...
    
    def __call__(self, values, clip=..., out=..., invalid=...): # -> NDArray[Unknown]:
        """
        Transform values using this stretch.

        Parameters
        ----------
        values : array-like
            The input values, which should already be normalized to the
            [0:1] range.
        clip : bool, optional
            If `True` (default), values outside the [0:1] range are
            clipped to the [0:1] range.
        out : ndarray, optional
            If specified, the output values will be placed in this array
            (typically used for in-place calculations).
        invalid : None or float, optional
            Value to assign NaN values generated by this class.  NaNs in
            the input ``values`` array are not changed.  This option is
            generally used with matplotlib normalization classes, where
            the ``invalid`` value should map to the matplotlib colormap
            "under" value (i.e., any finite value < 0).  If `None`, then
            NaN values are not replaced.  This keyword has no effect if
            ``clip=True``.

        Returns
        -------
        result : ndarray
            The transformed values.
        """
        ...
    
    @property
    def inverse(self): # -> PowerStretch:
        """A stretch object that performs the inverse operation."""
        ...
    


class PowerDistStretch(BaseStretch):
    r"""
    An alternative power stretch.

    The stretch is given by:

    .. math::
        y = \frac{a^x - 1}{a - 1}

    Parameters
    ----------
    a : float, optional
        The ``a`` parameter used in the above formula.  ``a`` must be
        greater than or equal to 0, but cannot be set to 1.  Default is
        1000.
    """
    def __init__(self, a=...) -> None:
        ...
    
    def __call__(self, values, clip=..., out=...): # -> NDArray[Unknown]:
        ...
    
    @property
    def inverse(self): # -> InvertedPowerDistStretch:
        """A stretch object that performs the inverse operation."""
        ...
    


class InvertedPowerDistStretch(BaseStretch):
    r"""
    Inverse transformation for
    `~astropy.image.scaling.PowerDistStretch`.

    The stretch is given by:

    .. math::
        y = \frac{\log(y (a-1) + 1)}{\log a}

    Parameters
    ----------
    a : float, optional
        The ``a`` parameter used in the above formula.  ``a`` must be
        greater than or equal to 0, but cannot be set to 1.  Default is
        1000.
    """
    def __init__(self, a=...) -> None:
        ...
    
    def __call__(self, values, clip=..., out=...): # -> NDArray[Unknown]:
        ...
    
    @property
    def inverse(self): # -> PowerDistStretch:
        """A stretch object that performs the inverse operation."""
        ...
    


class SquaredStretch(PowerStretch):
    r"""
    A convenience class for a power stretch of 2.

    The stretch is given by:

    .. math::
        y = x^2
    """
    def __init__(self) -> None:
        ...
    
    @property
    def inverse(self): # -> SqrtStretch:
        """A stretch object that performs the inverse operation."""
        ...
    


class LogStretch(BaseStretch):
    r"""
    A log stretch.

    The stretch is given by:

    .. math::
        y = \frac{\log{(a x + 1)}}{\log{(a + 1)}}

    Parameters
    ----------
    a : float
        The ``a`` parameter used in the above formula.  ``a`` must be
        greater than 0.  Default is 1000.
    """
    def __init__(self, a=...) -> None:
        ...
    
    def __call__(self, values, clip=..., out=..., invalid=...): # -> NDArray[Unknown]:
        """
        Transform values using this stretch.

        Parameters
        ----------
        values : array-like
            The input values, which should already be normalized to the
            [0:1] range.
        clip : bool, optional
            If `True` (default), values outside the [0:1] range are
            clipped to the [0:1] range.
        out : ndarray, optional
            If specified, the output values will be placed in this array
            (typically used for in-place calculations).
        invalid : None or float, optional
            Value to assign NaN values generated by this class.  NaNs in
            the input ``values`` array are not changed.  This option is
            generally used with matplotlib normalization classes, where
            the ``invalid`` value should map to the matplotlib colormap
            "under" value (i.e., any finite value < 0).  If `None`, then
            NaN values are not replaced.  This keyword has no effect if
            ``clip=True``.

        Returns
        -------
        result : ndarray
            The transformed values.
        """
        ...
    
    @property
    def inverse(self): # -> InvertedLogStretch:
        """A stretch object that performs the inverse operation."""
        ...
    


class InvertedLogStretch(BaseStretch):
    r"""
    Inverse transformation for `~astropy.image.scaling.LogStretch`.

    The stretch is given by:

    .. math::
        y = \frac{e^{y \log{a + 1}} - 1}{a} \\
        y = \frac{e^{y} (a + 1) - 1}{a}

    Parameters
    ----------
    a : float, optional
        The ``a`` parameter used in the above formula.  ``a`` must be
        greater than 0.  Default is 1000.
    """
    def __init__(self, a) -> None:
        ...
    
    def __call__(self, values, clip=..., out=...): # -> NDArray[Unknown]:
        ...
    
    @property
    def inverse(self): # -> LogStretch:
        """A stretch object that performs the inverse operation."""
        ...
    


class AsinhStretch(BaseStretch):
    r"""
    An asinh stretch.

    The stretch is given by:

    .. math::
        y = \frac{{\rm asinh}(x / a)}{{\rm asinh}(1 / a)}.

    Parameters
    ----------
    a : float, optional
        The ``a`` parameter used in the above formula.  The value of
        this parameter is where the asinh curve transitions from linear
        to logarithmic behavior, expressed as a fraction of the
        normalized image.  ``a`` must be greater than 0 and less than or
        equal to 1 (0 < a <= 1).  Default is 0.1.
    """
    def __init__(self, a=...) -> None:
        ...
    
    def __call__(self, values, clip=..., out=...): # -> NDArray[Unknown]:
        ...
    
    @property
    def inverse(self): # -> SinhStretch:
        """A stretch object that performs the inverse operation."""
        ...
    


class SinhStretch(BaseStretch):
    r"""
    A sinh stretch.

    The stretch is given by:

    .. math::
        y = \frac{{\rm sinh}(x / a)}{{\rm sinh}(1 / a)}

    Parameters
    ----------
    a : float, optional
        The ``a`` parameter used in the above formula.  ``a`` must be
        greater than 0 and less than or equal to 1 (0 < a <= 1).
        Default is 1/3.
    """
    def __init__(self, a=...) -> None:
        ...
    
    def __call__(self, values, clip=..., out=...): # -> NDArray[Unknown]:
        ...
    
    @property
    def inverse(self): # -> AsinhStretch:
        """A stretch object that performs the inverse operation."""
        ...
    


class HistEqStretch(BaseStretch):
    """
    A histogram equalization stretch.

    Parameters
    ----------
    data : array-like
        The data defining the equalization.
    values : array-like, optional
        The input image values, which should already be normalized to
        the [0:1] range.
    """
    def __init__(self, data, values=...) -> None:
        ...
    
    def __call__(self, values, clip=..., out=...): # -> NDArray[Unknown]:
        ...
    
    @property
    def inverse(self): # -> InvertedHistEqStretch:
        """A stretch object that performs the inverse operation."""
        ...
    


class InvertedHistEqStretch(BaseStretch):
    """
    Inverse transformation for `~astropy.image.scaling.HistEqStretch`.

    Parameters
    ----------
    data : array-like
        The data defining the equalization.
    values : array-like, optional
        The input image values, which should already be normalized to
        the [0:1] range.
    """
    def __init__(self, data, values=...) -> None:
        ...
    
    def __call__(self, values, clip=..., out=...): # -> NDArray[Unknown]:
        ...
    
    @property
    def inverse(self): # -> HistEqStretch:
        """A stretch object that performs the inverse operation."""
        ...
    


class ContrastBiasStretch(BaseStretch):
    r"""
    A stretch that takes into account contrast and bias.

    The stretch is given by:

    .. math::
        y = (x - {\rm bias}) * {\rm contrast} + 0.5

    and the output values are clipped to the [0:1] range.

    Parameters
    ----------
    contrast : float
        The contrast parameter (see the above formula).

    bias : float
        The bias parameter (see the above formula).
    """
    def __init__(self, contrast, bias) -> None:
        ...
    
    def __call__(self, values, clip=..., out=...): # -> NDArray[Unknown]:
        ...
    
    @property
    def inverse(self): # -> InvertedContrastBiasStretch:
        """A stretch object that performs the inverse operation."""
        ...
    


class InvertedContrastBiasStretch(BaseStretch):
    """
    Inverse transformation for ContrastBiasStretch.

    Parameters
    ----------
    contrast : float
        The contrast parameter (see
        `~astropy.visualization.ConstrastBiasStretch).

    bias : float
        The bias parameter (see
        `~astropy.visualization.ConstrastBiasStretch).
    """
    def __init__(self, contrast, bias) -> None:
        ...
    
    def __call__(self, values, clip=..., out=...): # -> NDArray[Unknown]:
        ...
    
    @property
    def inverse(self): # -> ContrastBiasStretch:
        """A stretch object that performs the inverse operation."""
        ...
    


class CompositeStretch(CompositeTransform, BaseStretch):
    """
    A combination of two stretches.

    Parameters
    ----------
    stretch_1 : :class:`astropy.visualization.BaseStretch`
        The first stretch to apply.
    stretch_2 : :class:`astropy.visualization.BaseStretch`
        The second stretch to apply.
    """
    def __call__(self, values, clip=..., out=...):
        ...
    


