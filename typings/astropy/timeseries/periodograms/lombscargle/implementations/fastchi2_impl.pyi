"""
This type stub file was generated by pyright.
"""

def lombscargle_fastchi2(t, y, dy, f0, df, Nf, normalization=..., fit_mean=..., center_data=..., nterms=..., use_fft=..., trig_sum_kwds=...): # -> NDArray[floating[Any]] | Any:
    """Lomb-Scargle Periodogram

    This implements a fast chi-squared periodogram using the algorithm
    outlined in [4]_. The result is identical to the standard Lomb-Scargle
    periodogram. The advantage of this algorithm is the
    ability to compute multiterm periodograms relatively quickly.

    Parameters
    ----------
    t, y, dy : array-like
        times, values, and errors of the data points. These should be
        broadcastable to the same shape. None should be `~astropy.units.Quantity`.
    f0, df, Nf : (float, float, int)
        parameters describing the frequency grid, f = f0 + df * arange(Nf).
    normalization : str, optional
        Normalization to use for the periodogram.
        Options are 'standard', 'model', 'log', or 'psd'.
    fit_mean : bool, optional
        if True, include a constant offset as part of the model at each
        frequency. This can lead to more accurate results, especially in the
        case of incomplete phase coverage.
    center_data : bool, optional
        if True, pre-center the data by subtracting the weighted mean
        of the input data. This is especially important if ``fit_mean = False``
    nterms : int, optional
        Number of Fourier terms in the fit

    Returns
    -------
    power : array-like
        Lomb-Scargle power associated with each frequency.
        Units of the result depend on the normalization.

    References
    ----------
    .. [1] M. Zechmeister and M. Kurster, A&A 496, 577-584 (2009)
    .. [2] W. Press et al, Numerical Recipes in C (2002)
    .. [3] Scargle, J.D. ApJ 263:835-853 (1982)
    .. [4] Palmer, J. ApJ 695:496-502 (2009)
    """
    ...

