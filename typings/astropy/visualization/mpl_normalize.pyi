"""
This type stub file was generated by pyright.
"""

"""
Normalization class for Matplotlib that can be used to produce
colorbars.
"""
__all__ = ['ImageNormalize', 'simple_norm', 'imshow_norm']
__doctest_requires__ = ...
class ImageNormalize(Normalize):
    """
    Normalization class to be used with Matplotlib.

    Parameters
    ----------
    data : ndarray, optional
        The image array.  This input is used only if ``interval`` is
        also input.  ``data`` and ``interval`` are used to compute the
        vmin and/or vmax values only if ``vmin`` or ``vmax`` are not
        input.
    interval : `~astropy.visualization.BaseInterval` subclass instance, optional
        The interval object to apply to the input ``data`` to determine
        the ``vmin`` and ``vmax`` values.  This input is used only if
        ``data`` is also input.  ``data`` and ``interval`` are used to
        compute the vmin and/or vmax values only if ``vmin`` or ``vmax``
        are not input.
    vmin, vmax : float, optional
        The minimum and maximum levels to show for the data.  The
        ``vmin`` and ``vmax`` inputs override any calculated values from
        the ``interval`` and ``data`` inputs.
    stretch : `~astropy.visualization.BaseStretch` subclass instance
        The stretch object to apply to the data.  The default is
        `~astropy.visualization.LinearStretch`.
    clip : bool, optional
        If `True`, data values outside the [0:1] range are clipped to
        the [0:1] range.
    invalid : None or float, optional
        Value to assign NaN values generated by this class.  NaNs in the
        input ``data`` array are not changed.  For matplotlib
        normalization, the ``invalid`` value should map to the
        matplotlib colormap "under" value (i.e., any finite value < 0).
        If `None`, then NaN values are not replaced.  This keyword has
        no effect if ``clip=True``.
    """
    def __init__(self, data=..., interval=..., vmin=..., vmax=..., stretch=..., clip=..., invalid=...) -> None:
        ...
    
    def __call__(self, values, clip=..., invalid=...):
        """
        Transform values using this normalization.

        Parameters
        ----------
        values : array-like
            The input values.
        clip : bool, optional
            If `True`, values outside the [0:1] range are clipped to the
            [0:1] range.  If `None` then the ``clip`` value from the
            `ImageNormalize` instance is used (the default of which is
            `False`).
        invalid : None or float, optional
            Value to assign NaN values generated by this class.  NaNs in
            the input ``data`` array are not changed.  For matplotlib
            normalization, the ``invalid`` value should map to the
            matplotlib colormap "under" value (i.e., any finite value <
            0).  If `None`, then the `ImageNormalize` instance value is
            used.  This keyword has no effect if ``clip=True``.
        """
        ...
    
    def inverse(self, values, invalid=...):
        ...
    


def simple_norm(data, stretch=..., power=..., asinh_a=..., min_cut=..., max_cut=..., min_percent=..., max_percent=..., percent=..., clip=..., log_a=..., invalid=...):
    """
    Return a Normalization class that can be used for displaying images
    with Matplotlib.

    This function enables only a subset of image stretching functions
    available in `~astropy.visualization.mpl_normalize.ImageNormalize`.

    This function is used by the
    ``astropy.visualization.scripts.fits2bitmap`` script.

    Parameters
    ----------
    data : ndarray
        The image array.

    stretch : {'linear', 'sqrt', 'power', log', 'asinh'}, optional
        The stretch function to apply to the image.  The default is
        'linear'.

    power : float, optional
        The power index for ``stretch='power'``.  The default is 1.0.

    asinh_a : float, optional
        For ``stretch='asinh'``, the value where the asinh curve
        transitions from linear to logarithmic behavior, expressed as a
        fraction of the normalized image.  Must be in the range between
        0 and 1.  The default is 0.1.

    min_cut : float, optional
        The pixel value of the minimum cut level.  Data values less than
        ``min_cut`` will set to ``min_cut`` before stretching the image.
        The default is the image minimum.  ``min_cut`` overrides
        ``min_percent``.

    max_cut : float, optional
        The pixel value of the maximum cut level.  Data values greater
        than ``min_cut`` will set to ``min_cut`` before stretching the
        image.  The default is the image maximum.  ``max_cut`` overrides
        ``max_percent``.

    min_percent : float, optional
        The percentile value used to determine the pixel value of
        minimum cut level.  The default is 0.0.  ``min_percent``
        overrides ``percent``.

    max_percent : float, optional
        The percentile value used to determine the pixel value of
        maximum cut level.  The default is 100.0.  ``max_percent``
        overrides ``percent``.

    percent : float, optional
        The percentage of the image values used to determine the pixel
        values of the minimum and maximum cut levels.  The lower cut
        level will set at the ``(100 - percent) / 2`` percentile, while
        the upper cut level will be set at the ``(100 + percent) / 2``
        percentile.  The default is 100.0.  ``percent`` is ignored if
        either ``min_percent`` or ``max_percent`` is input.

    clip : bool, optional
        If `True`, data values outside the [0:1] range are clipped to
        the [0:1] range.

    log_a : float, optional
        The log index for ``stretch='log'``. The default is 1000.

    invalid : None or float, optional
        Value to assign NaN values generated by the normalization.  NaNs
        in the input ``data`` array are not changed.  For matplotlib
        normalization, the ``invalid`` value should map to the
        matplotlib colormap "under" value (i.e., any finite value < 0).
        If `None`, then NaN values are not replaced.  This keyword has
        no effect if ``clip=True``.

    Returns
    -------
    result : `ImageNormalize` instance
        An `ImageNormalize` instance that can be used for displaying
        images with Matplotlib.
    """
    ...

_norm_sig = ...
def imshow_norm(data, ax=..., **kwargs):
    """ A convenience function to call matplotlib's `matplotlib.pyplot.imshow`
    function, using an `ImageNormalize` object as the normalization.

    Parameters
    ----------
    data : 2D or 3D array-like
        The data to show. Can be whatever `~matplotlib.pyplot.imshow` and
        `ImageNormalize` both accept. See `~matplotlib.pyplot.imshow`.
    ax : None or `~matplotlib.axes.Axes`, optional
        If None, use pyplot's imshow.  Otherwise, calls ``imshow`` method of
        the supplied axes.
    **kwargs : dict, optional
        All other keyword arguments are parsed first by the
        `ImageNormalize` initializer, then to
        `~matplotlib.pyplot.imshow`.

    Returns
    -------
    result : tuple
        A tuple containing the `~matplotlib.image.AxesImage` generated
        by `~matplotlib.pyplot.imshow` as well as the `ImageNormalize`
        instance.

    Notes
    -----
    The ``norm`` matplotlib keyword is not supported.

    Examples
    --------
    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt
        from astropy.visualization import (imshow_norm, MinMaxInterval,
                                           SqrtStretch)

        # Generate and display a test image
        image = np.arange(65536).reshape((256, 256))
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        im, norm = imshow_norm(image, ax, origin='lower',
                               interval=MinMaxInterval(),
                               stretch=SqrtStretch())
        fig.colorbar(im)
    """
    ...
