"""
This type stub file was generated by pyright.
"""

from .core import Model

"""
Tabular models.

Tabular models of any dimension can be created using `tabular_model`.
For convenience `Tabular1D` and `Tabular2D` are provided.

Examples
--------
>>> table = np.array([[ 3.,  0.,  0.],
...                  [ 0.,  2.,  0.],
...                  [ 0.,  0.,  0.]])
>>> points = ([1, 2, 3], [1, 2, 3])
>>> t2 = Tabular2D(points, lookup_table=table, bounds_error=False,
...                fill_value=None, method='nearest')

"""
__all__ = ['tabular_model', 'Tabular1D', 'Tabular2D']
__doctest_requires__ = ...
class _Tabular(Model):
    """
    Returns an interpolated lookup table value.

    Parameters
    ----------
    points : tuple of ndarray of float, optional
        The points defining the regular grid in n dimensions.
        ndarray must have shapes (m1, ), ..., (mn, ),
    lookup_table : array-like
        The data on a regular grid in n dimensions.
        Must have shapes (m1, ..., mn, ...)
    method : str, optional
        The method of interpolation to perform. Supported are "linear" and
        "nearest", and "splinef2d". "splinef2d" is only supported for
        2-dimensional data. Default is "linear".
    bounds_error : bool, optional
        If True, when interpolated values are requested outside of the
        domain of the input data, a ValueError is raised.
        If False, then ``fill_value`` is used.
    fill_value : float or `~astropy.units.Quantity`, optional
        If provided, the value to use for points outside of the
        interpolation domain. If None, values outside
        the domain are extrapolated.  Extrapolation is not supported by method
        "splinef2d". If Quantity is given, it will be converted to the unit of
        ``lookup_table``, if applicable.

    Returns
    -------
    value : ndarray
        Interpolated values at input coordinates.

    Raises
    ------
    ImportError
        Scipy is not installed.

    Notes
    -----
    Uses `scipy.interpolate.interpn`.

    """
    linear = ...
    fittable = ...
    standard_broadcasting = ...
    _is_dynamic = ...
    _id = ...
    def __init__(self, points=..., lookup_table=..., method=..., bounds_error=..., fill_value=..., **kwargs) -> None:
        ...
    
    def __repr__(self): # -> str:
        ...
    
    def __str__(self) -> str:
        ...
    
    @property
    def input_units(self): # -> dict[str | Unknown, <subclass of Unit and StructuredUnit> | StructuredUnit | Unit | UnitBase | Unknown | None] | None:
        ...
    
    @property
    def return_units(self): # -> dict[str | Any | Unknown, <subclass of Unit and StructuredUnit> | StructuredUnit | Unit | UnitBase | None] | None:
        ...
    
    @property
    def bounding_box(self): # -> tuple[Any, Any] | list[tuple[Any, Any]]:
        """
        Tuple defining the default ``bounding_box`` limits,
        ``(points_low, points_high)``.

        Examples
        --------
        >>> from astropy.modeling.models import Tabular1D, Tabular2D
        >>> t1 = Tabular1D(points=[1, 2, 3], lookup_table=[10, 20, 30])
        >>> t1.bounding_box
        ModelBoundingBox(
            intervals={
                x: Interval(lower=1, upper=3)
            }
            model=Tabular1D(inputs=('x',))
            order='C'
        )
        >>> t2 = Tabular2D(points=[[1, 2, 3], [2, 3, 4]],
        ...                lookup_table=[[10, 20, 30], [20, 30, 40]])
        >>> t2.bounding_box
        ModelBoundingBox(
            intervals={
                x: Interval(lower=1, upper=3)
                y: Interval(lower=2, upper=4)
            }
            model=Tabular2D(inputs=('x', 'y'))
            order='C'
        )

        """
        ...
    
    def evaluate(self, *inputs): # -> list[Unknown]:
        """
        Return the interpolated values at the input coordinates.

        Parameters
        ----------
        inputs : list of scalar or list of ndarray
            Input coordinates. The number of inputs must be equal
            to the dimensions of the lookup table.
        """
        ...
    
    @property
    def inverse(self): # -> Tabular1D:
        ...
    


def tabular_model(dim, name=...): # -> Type[_]:
    """
    Make a ``Tabular`` model where ``n_inputs`` is
    based on the dimension of the lookup_table.

    This model has to be further initialized and when evaluated
    returns the interpolated values.

    Parameters
    ----------
    dim : int
        Dimensions of the lookup table.
    name : str
        Name for the class.

    Examples
    --------
    >>> table = np.array([[3., 0., 0.],
    ...                   [0., 2., 0.],
    ...                   [0., 0., 0.]])

    >>> tab = tabular_model(2, name='Tabular2D')
    >>> print(tab)
    <class 'astropy.modeling.tabular.Tabular2D'>
    Name: Tabular2D
    N_inputs: 2
    N_outputs: 1

    >>> points = ([1, 2, 3], [1, 2, 3])

    Setting fill_value to None, allows extrapolation.
    >>> m = tab(points, lookup_table=table, name='my_table',
    ...         bounds_error=False, fill_value=None, method='nearest')

    >>> xinterp = [0, 1, 1.5, 2.72, 3.14]
    >>> m(xinterp, xinterp)  # doctest: +FLOAT_CMP
    array([3., 3., 3., 0., 0.])

    """
    ...

Tabular1D = ...
Tabular2D = ...
_tab_docs = ...
