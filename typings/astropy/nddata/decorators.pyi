"""
This type stub file was generated by pyright.
"""

__all__ = ['support_nddata']
SUPPORTED_PROPERTIES = ...
def support_nddata(_func=..., accepts=..., repack=..., returns=..., keeps=..., **attribute_argument_mapping): # -> ((data: Unknown, *args: Unknown, **kwargs: Unknown) -> Unknown) | ((func: Unknown) -> ((data: Unknown, *args: Unknown, **kwargs: Unknown) -> Unknown)):
    """Decorator to wrap functions that could accept an NDData instance with
    its properties passed as function arguments.

    Parameters
    ----------
    _func : callable, None, optional
        The function to decorate or ``None`` if used as factory. The first
        positional argument should be ``data`` and take a numpy array. It is
        possible to overwrite the name, see ``attribute_argument_mapping``
        argument.
        Default is ``None``.

    accepts : class, optional
        The class or subclass of ``NDData`` that should be unpacked before
        calling the function.
        Default is ``NDData``

    repack : bool, optional
        Should be ``True`` if the return should be converted to the input
        class again after the wrapped function call.
        Default is ``False``.

        .. note::
           Must be ``True`` if either one of ``returns`` or ``keeps``
           is specified.

    returns : iterable, None, optional
        An iterable containing strings which returned value should be set
        on the class. For example if a function returns data and mask, this
        should be ``['data', 'mask']``. If ``None`` assume the function only
        returns one argument: ``'data'``.
        Default is ``None``.

        .. note::
           Must be ``None`` if ``repack=False``.

    keeps : iterable. None, optional
        An iterable containing strings that indicate which values should be
        copied from the original input to the returned class. If ``None``
        assume that no attributes are copied.
        Default is ``None``.

        .. note::
           Must be ``None`` if ``repack=False``.

    attribute_argument_mapping :
        Keyword parameters that optionally indicate which function argument
        should be interpreted as which attribute on the input. By default
        it assumes the function takes a ``data`` argument as first argument,
        but if the first argument is called ``input`` one should pass
        ``support_nddata(..., data='input')`` to the function.

    Returns
    -------
    decorator_factory or decorated_function : callable
        If ``_func=None`` this returns a decorator, otherwise it returns the
        decorated ``_func``.

    Notes
    -----
    If properties of ``NDData`` are set but have no corresponding function
    argument a Warning is shown.

    If a property is set of the ``NDData`` are set and an explicit argument is
    given, the explicitly given argument is used and a Warning is shown.

    The supported properties are:

    - ``mask``
    - ``unit``
    - ``wcs``
    - ``meta``
    - ``uncertainty``
    - ``flags``

    Examples
    --------

    This function takes a Numpy array for the data, and some WCS information
    with the ``wcs`` keyword argument::

        def downsample(data, wcs=None):
            # downsample data and optionally WCS here
            pass

    However, you might have an NDData instance that has the ``wcs`` property
    set and you would like to be able to call the function with
    ``downsample(my_nddata)`` and have the WCS information, if present,
    automatically be passed to the ``wcs`` keyword argument.

    This decorator can be used to make this possible::

        @support_nddata
        def downsample(data, wcs=None):
            # downsample data and optionally WCS here
            pass

    This function can now either be called as before, specifying the data and
    WCS separately, or an NDData instance can be passed to the ``data``
    argument.
    """
    ...

