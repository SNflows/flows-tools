"""
This type stub file was generated by pyright.
"""

"""
Main Lomb-Scargle Implementation

The ``lombscargle`` function here is essentially a sophisticated switch
statement for the various implementations available in this submodule
"""
__all__ = ['lombscargle', 'available_methods']
METHODS = ...
def available_methods(): # -> list[str]:
    ...

def validate_method(method, dy, fit_mean, nterms, frequency, assume_regular_frequency): # -> Literal['fastchi2', 'chi2', 'fast', 'scipy', 'cython']:
    """
    Validate the method argument, and if method='auto'
    choose the appropriate method
    """
    ...

def lombscargle(t, y, dy=..., frequency=..., method=..., assume_regular_frequency=..., normalization=..., fit_mean=..., center_data=..., method_kwds=..., nterms=...): # -> Any | ndarray[Any, dtype[floating[Any]]]:
    """
    Compute the Lomb-scargle Periodogram with a given method.

    Parameters
    ----------
    t : array-like
        sequence of observation times
    y : array-like
        sequence of observations associated with times t
    dy : float or array-like, optional
        error or sequence of observational errors associated with times t
    frequency : array-like
        frequencies (not angular frequencies) at which to evaluate the
        periodogram. If not specified, optimal frequencies will be chosen using
        a heuristic which will attempt to provide sufficient frequency range
        and sampling so that peaks will not be missed. Note that in order to
        use method='fast', frequencies must be regularly spaced.
    method : str, optional
        specify the lomb scargle implementation to use. Options are:

        - 'auto': choose the best method based on the input
        - 'fast': use the O[N log N] fast method. Note that this requires
          evenly-spaced frequencies: by default this will be checked unless
          ``assume_regular_frequency`` is set to True.
        - `slow`: use the O[N^2] pure-python implementation
        - `chi2`: use the O[N^2] chi2/linear-fitting implementation
        - `fastchi2`: use the O[N log N] chi2 implementation. Note that this
          requires evenly-spaced frequencies: by default this will be checked
          unless `assume_regular_frequency` is set to True.
        - `scipy`: use ``scipy.signal.lombscargle``, which is an O[N^2]
          implementation written in C. Note that this does not support
          heteroskedastic errors.

    assume_regular_frequency : bool, optional
        if True, assume that the input frequency is of the form
        freq = f0 + df * np.arange(N). Only referenced if method is 'auto'
        or 'fast'.
    normalization : str, optional
        Normalization to use for the periodogram.
        Options are 'standard' or 'psd'.
    fit_mean : bool, optional
        if True, include a constant offset as part of the model at each
        frequency. This can lead to more accurate results, especially in the
        case of incomplete phase coverage.
    center_data : bool, optional
        if True, pre-center the data by subtracting the weighted mean
        of the input data. This is especially important if `fit_mean = False`
    method_kwds : dict, optional
        additional keywords to pass to the lomb-scargle method
    nterms : int, optional
        number of Fourier terms to use in the periodogram.
        Not supported with every method.

    Returns
    -------
    PLS : array-like
        Lomb-Scargle power associated with each frequency omega
    """
    ...
