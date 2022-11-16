"""
This type stub file was generated by pyright.
"""

__doctest_skip__ = ...
__all__ = ['obsgeo_to_frame', 'add_stokes_axis_to_wcs', 'celestial_frame_to_wcs', 'wcs_to_celestial_frame', 'proj_plane_pixel_scales', 'proj_plane_pixel_area', 'is_proj_plane_distorted', 'non_celestial_pixel_scales', 'skycoord_to_pixel', 'pixel_to_skycoord', 'custom_wcs_to_frame_mappings', 'custom_frame_to_wcs_mappings', 'pixel_to_pixel', 'local_partial_pixel_derivatives', 'fit_wcs_from_points']
def add_stokes_axis_to_wcs(wcs, add_before_ind):
    """
    Add a new Stokes axis that is uncorrelated with any other axes.

    Parameters
    ----------
    wcs : `~astropy.wcs.WCS`
        The WCS to add to
    add_before_ind : int
        Index of the WCS to insert the new Stokes axis in front of.
        To add at the end, do add_before_ind = wcs.wcs.naxis
        The beginning is at position 0.

    Returns
    -------
    `~astropy.wcs.WCS`
        A new `~astropy.wcs.WCS` instance with an additional axis
    """
    ...

WCS_FRAME_MAPPINGS = ...
FRAME_WCS_MAPPINGS = ...
class custom_wcs_to_frame_mappings:
    def __init__(self, mappings=...) -> None:
        ...
    
    def __enter__(self): # -> None:
        ...
    
    def __exit__(self, type, value, tb): # -> None:
        ...
    


custom_frame_mappings = custom_wcs_to_frame_mappings
class custom_frame_to_wcs_mappings:
    def __init__(self, mappings=...) -> None:
        ...
    
    def __enter__(self): # -> None:
        ...
    
    def __exit__(self, type, value, tb): # -> None:
        ...
    


def wcs_to_celestial_frame(wcs): # -> FK4 | FK4NoETerms | FK5 | ICRS | Galactic | ITRS:
    """
    For a given WCS, return the coordinate frame that matches the celestial
    component of the WCS.

    Parameters
    ----------
    wcs : :class:`~astropy.wcs.WCS` instance
        The WCS to find the frame for

    Returns
    -------
    frame : :class:`~astropy.coordinates.baseframe.BaseCoordinateFrame` subclass instance
        An instance of a :class:`~astropy.coordinates.baseframe.BaseCoordinateFrame`
        subclass instance that best matches the specified WCS.

    Notes
    -----

    To extend this function to frames not defined in astropy.coordinates, you
    can write your own function which should take a :class:`~astropy.wcs.WCS`
    instance and should return either an instance of a frame, or `None` if no
    matching frame was found. You can register this function temporarily with::

        >>> from astropy.wcs.utils import wcs_to_celestial_frame, custom_wcs_to_frame_mappings
        >>> with custom_wcs_to_frame_mappings(my_function):
        ...     wcs_to_celestial_frame(...)

    """
    ...

def celestial_frame_to_wcs(frame, projection=...): # -> WCS:
    """
    For a given coordinate frame, return the corresponding WCS object.

    Note that the returned WCS object has only the elements corresponding to
    coordinate frames set (e.g. ctype, equinox, radesys).

    Parameters
    ----------
    frame : :class:`~astropy.coordinates.baseframe.BaseCoordinateFrame` subclass instance
        An instance of a :class:`~astropy.coordinates.baseframe.BaseCoordinateFrame`
        subclass instance for which to find the WCS
    projection : str
        Projection code to use in ctype, if applicable

    Returns
    -------
    wcs : :class:`~astropy.wcs.WCS` instance
        The corresponding WCS object

    Examples
    --------

    ::

        >>> from astropy.wcs.utils import celestial_frame_to_wcs
        >>> from astropy.coordinates import FK5
        >>> frame = FK5(equinox='J2010')
        >>> wcs = celestial_frame_to_wcs(frame)
        >>> wcs.to_header()
        WCSAXES =                    2 / Number of coordinate axes
        CRPIX1  =                  0.0 / Pixel coordinate of reference point
        CRPIX2  =                  0.0 / Pixel coordinate of reference point
        CDELT1  =                  1.0 / [deg] Coordinate increment at reference point
        CDELT2  =                  1.0 / [deg] Coordinate increment at reference point
        CUNIT1  = 'deg'                / Units of coordinate increment and value
        CUNIT2  = 'deg'                / Units of coordinate increment and value
        CTYPE1  = 'RA---TAN'           / Right ascension, gnomonic projection
        CTYPE2  = 'DEC--TAN'           / Declination, gnomonic projection
        CRVAL1  =                  0.0 / [deg] Coordinate value at reference point
        CRVAL2  =                  0.0 / [deg] Coordinate value at reference point
        LONPOLE =                180.0 / [deg] Native longitude of celestial pole
        LATPOLE =                  0.0 / [deg] Native latitude of celestial pole
        RADESYS = 'FK5'                / Equatorial coordinate system
        EQUINOX =               2010.0 / [yr] Equinox of equatorial coordinates


    Notes
    -----

    To extend this function to frames not defined in astropy.coordinates, you
    can write your own function which should take a
    :class:`~astropy.coordinates.baseframe.BaseCoordinateFrame` subclass
    instance and a projection (given as a string) and should return either a WCS
    instance, or `None` if the WCS could not be determined. You can register
    this function temporarily with::

        >>> from astropy.wcs.utils import celestial_frame_to_wcs, custom_frame_to_wcs_mappings
        >>> with custom_frame_to_wcs_mappings(my_function):
        ...     celestial_frame_to_wcs(...)

    """
    ...

def proj_plane_pixel_scales(wcs): # -> Any:
    """
    For a WCS returns pixel scales along each axis of the image pixel at
    the ``CRPIX`` location once it is projected onto the
    "plane of intermediate world coordinates" as defined in
    `Greisen & Calabretta 2002, A&A, 395, 1061 <https://ui.adsabs.harvard.edu/abs/2002A%26A...395.1061G>`_.

    .. note::
        This function is concerned **only** about the transformation
        "image plane"->"projection plane" and **not** about the
        transformation "celestial sphere"->"projection plane"->"image plane".
        Therefore, this function ignores distortions arising due to
        non-linear nature of most projections.

    .. note::
        In order to compute the scales corresponding to celestial axes only,
        make sure that the input `~astropy.wcs.WCS` object contains
        celestial axes only, e.g., by passing in the
        `~astropy.wcs.WCS.celestial` WCS object.

    Parameters
    ----------
    wcs : `~astropy.wcs.WCS`
        A world coordinate system object.

    Returns
    -------
    scale : ndarray
        A vector (`~numpy.ndarray`) of projection plane increments
        corresponding to each pixel side (axis). The units of the returned
        results are the same as the units of `~astropy.wcs.Wcsprm.cdelt`,
        `~astropy.wcs.Wcsprm.crval`, and `~astropy.wcs.Wcsprm.cd` for
        the celestial WCS and can be obtained by inquiring the value
        of `~astropy.wcs.Wcsprm.cunit` property of the input
        `~astropy.wcs.WCS` WCS object.

    See Also
    --------
    astropy.wcs.utils.proj_plane_pixel_area

    """
    ...

def proj_plane_pixel_area(wcs): # -> Any:
    """
    For a **celestial** WCS (see `astropy.wcs.WCS.celestial`) returns pixel
    area of the image pixel at the ``CRPIX`` location once it is projected
    onto the "plane of intermediate world coordinates" as defined in
    `Greisen & Calabretta 2002, A&A, 395, 1061 <https://ui.adsabs.harvard.edu/abs/2002A%26A...395.1061G>`_.

    .. note::
        This function is concerned **only** about the transformation
        "image plane"->"projection plane" and **not** about the
        transformation "celestial sphere"->"projection plane"->"image plane".
        Therefore, this function ignores distortions arising due to
        non-linear nature of most projections.

    .. note::
        In order to compute the area of pixels corresponding to celestial
        axes only, this function uses the `~astropy.wcs.WCS.celestial` WCS
        object of the input ``wcs``.  This is different from the
        `~astropy.wcs.utils.proj_plane_pixel_scales` function
        that computes the scales for the axes of the input WCS itself.

    Parameters
    ----------
    wcs : `~astropy.wcs.WCS`
        A world coordinate system object.

    Returns
    -------
    area : float
        Area (in the projection plane) of the pixel at ``CRPIX`` location.
        The units of the returned result are the same as the units of
        the `~astropy.wcs.Wcsprm.cdelt`, `~astropy.wcs.Wcsprm.crval`,
        and `~astropy.wcs.Wcsprm.cd` for the celestial WCS and can be
        obtained by inquiring the value of `~astropy.wcs.Wcsprm.cunit`
        property of the `~astropy.wcs.WCS.celestial` WCS object.

    Raises
    ------
    ValueError
        Pixel area is defined only for 2D pixels. Most likely the
        `~astropy.wcs.Wcsprm.cd` matrix of the `~astropy.wcs.WCS.celestial`
        WCS is not a square matrix of second order.

    Notes
    -----

    Depending on the application, square root of the pixel area can be used to
    represent a single pixel scale of an equivalent square pixel
    whose area is equal to the area of a generally non-square pixel.

    See Also
    --------
    astropy.wcs.utils.proj_plane_pixel_scales

    """
    ...

def is_proj_plane_distorted(wcs, maxerr=...): # -> bool:
    r"""
    For a WCS returns `False` if square image (detector) pixels stay square
    when projected onto the "plane of intermediate world coordinates"
    as defined in
    `Greisen & Calabretta 2002, A&A, 395, 1061 <https://ui.adsabs.harvard.edu/abs/2002A%26A...395.1061G>`_.
    It will return `True` if transformation from image (detector) coordinates
    to the focal plane coordinates is non-orthogonal or if WCS contains
    non-linear (e.g., SIP) distortions.

    .. note::
        Since this function is concerned **only** about the transformation
        "image plane"->"focal plane" and **not** about the transformation
        "celestial sphere"->"focal plane"->"image plane",
        this function ignores distortions arising due to non-linear nature
        of most projections.

    Let's denote by *C* either the original or the reconstructed
    (from ``PC`` and ``CDELT``) CD matrix. `is_proj_plane_distorted`
    verifies that the transformation from image (detector) coordinates
    to the focal plane coordinates is orthogonal using the following
    check:

    .. math::
        \left \| \frac{C \cdot C^{\mathrm{T}}}
        {| det(C)|} - I \right \|_{\mathrm{max}} < \epsilon .

    Parameters
    ----------
    wcs : `~astropy.wcs.WCS`
        World coordinate system object

    maxerr : float, optional
        Accuracy to which the CD matrix, **normalized** such
        that :math:`|det(CD)|=1`, should be close to being an
        orthogonal matrix as described in the above equation
        (see :math:`\epsilon`).

    Returns
    -------
    distorted : bool
        Returns `True` if focal (projection) plane is distorted and `False`
        otherwise.

    """
    ...

def non_celestial_pixel_scales(inwcs): # -> Any:
    """
    Calculate the pixel scale along each axis of a non-celestial WCS,
    for example one with mixed spectral and spatial axes.

    Parameters
    ----------
    inwcs : `~astropy.wcs.WCS`
        The world coordinate system object.

    Returns
    -------
    scale : `numpy.ndarray`
        The pixel scale along each axis.
    """
    ...

def skycoord_to_pixel(coords, wcs, origin=..., mode=...): # -> tuple[Unknown, Unknown]:
    """
    Convert a set of SkyCoord coordinates into pixels.

    Parameters
    ----------
    coords : `~astropy.coordinates.SkyCoord`
        The coordinates to convert.
    wcs : `~astropy.wcs.WCS`
        The WCS transformation to use.
    origin : int
        Whether to return 0 or 1-based pixel coordinates.
    mode : 'all' or 'wcs'
        Whether to do the transformation including distortions (``'all'``) or
        only including only the core WCS transformation (``'wcs'``).

    Returns
    -------
    xp, yp : `numpy.ndarray`
        The pixel coordinates

    See Also
    --------
    astropy.coordinates.SkyCoord.from_pixel
    """
    ...

def pixel_to_skycoord(xp, yp, wcs, origin=..., mode=..., cls=...): # -> SkyCoord:
    """
    Convert a set of pixel coordinates into a `~astropy.coordinates.SkyCoord`
    coordinate.

    Parameters
    ----------
    xp, yp : float or ndarray
        The coordinates to convert.
    wcs : `~astropy.wcs.WCS`
        The WCS transformation to use.
    origin : int
        Whether to return 0 or 1-based pixel coordinates.
    mode : 'all' or 'wcs'
        Whether to do the transformation including distortions (``'all'``) or
        only including only the core WCS transformation (``'wcs'``).
    cls : class or None
        The class of object to create.  Should be a
        `~astropy.coordinates.SkyCoord` subclass.  If None, defaults to
        `~astropy.coordinates.SkyCoord`.

    Returns
    -------
    coords : `~astropy.coordinates.SkyCoord` subclass
        The celestial coordinates. Whatever ``cls`` type is.

    See Also
    --------
    astropy.coordinates.SkyCoord.from_pixel
    """
    ...

def pixel_to_pixel(wcs_in, wcs_out, *inputs):
    """
    Transform pixel coordinates in a dataset with a WCS to pixel coordinates
    in another dataset with a different WCS.

    This function is designed to efficiently deal with input pixel arrays that
    are broadcasted views of smaller arrays, and is compatible with any
    APE14-compliant WCS.

    Parameters
    ----------
    wcs_in : `~astropy.wcs.wcsapi.BaseHighLevelWCS`
        A WCS object for the original dataset which complies with the
        high-level shared APE 14 WCS API.
    wcs_out : `~astropy.wcs.wcsapi.BaseHighLevelWCS`
        A WCS object for the target dataset which complies with the
        high-level shared APE 14 WCS API.
    *inputs :
        Scalars or arrays giving the pixel coordinates to transform.
    """
    ...

def local_partial_pixel_derivatives(wcs, *pixel, normalize_by_world=...): # -> Any | NDArray[float64]:
    """
    Return a matrix of shape ``(world_n_dim, pixel_n_dim)`` where each entry
    ``[i, j]`` is the partial derivative d(world_i)/d(pixel_j) at the requested
    pixel position.

    Parameters
    ----------
    wcs : `~astropy.wcs.WCS`
        The WCS transformation to evaluate the derivatives for.
    *pixel : float
        The scalar pixel coordinates at which to evaluate the derivatives.
    normalize_by_world : bool
        If `True`, the matrix is normalized so that for each world entry
        the derivatives add up to 1.
    """
    ...

def fit_wcs_from_points(xy, world_coords, proj_point=..., projection=..., sip_degree=...): # -> WCS | str:
    """
    Given two matching sets of coordinates on detector and sky,
    compute the WCS.

    Fits a WCS object to matched set of input detector and sky coordinates.
    Optionally, a SIP can be fit to account for geometric
    distortion. Returns an `~astropy.wcs.WCS` object with the best fit
    parameters for mapping between input pixel and sky coordinates.

    The projection type (default 'TAN') can passed in as a string, one of
    the valid three-letter projection codes - or as a WCS object with
    projection keywords already set. Note that if an input WCS has any
    non-polynomial distortion, this will be applied and reflected in the
    fit terms and coefficients. Passing in a WCS object in this way essentially
    allows it to be refit based on the matched input coordinates and projection
    point, but take care when using this option as non-projection related
    keywords in the input might cause unexpected behavior.

    Notes
    -----
    - The fiducial point for the spherical projection can be set to 'center'
      to use the mean position of input sky coordinates, or as an
      `~astropy.coordinates.SkyCoord` object.
    - Units in all output WCS objects will always be in degrees.
    - If the coordinate frame differs between `~astropy.coordinates.SkyCoord`
      objects passed in for ``world_coords`` and ``proj_point``, the frame for
      ``world_coords``  will override as the frame for the output WCS.
    - If a WCS object is passed in to ``projection`` the CD/PC matrix will
      be used as an initial guess for the fit. If this is known to be
      significantly off and may throw off the fit, set to the identity matrix
      (for example, by doing wcs.wcs.pc = [(1., 0.,), (0., 1.)])

    Parameters
    ----------
    xy : (`numpy.ndarray`, `numpy.ndarray`) tuple
        x & y pixel coordinates.
    world_coords : `~astropy.coordinates.SkyCoord`
        Skycoord object with world coordinates.
    proj_point : 'center' or ~astropy.coordinates.SkyCoord`
        Defaults to 'center', in which the geometric center of input world
        coordinates will be used as the projection point. To specify an exact
        point for the projection, a Skycoord object with a coordinate pair can
        be passed in. For consistency, the units and frame of these coordinates
        will be transformed to match ``world_coords`` if they don't.
    projection : str or `~astropy.wcs.WCS`
        Three letter projection code, of any of standard projections defined
        in the FITS WCS standard. Optionally, a WCS object with projection
        keywords set may be passed in.
    sip_degree : None or int
        If set to a non-zero integer value, will fit SIP of degree
        ``sip_degree`` to model geometric distortion. Defaults to None, meaning
        no distortion corrections will be fit.

    Returns
    -------
    wcs : `~astropy.wcs.WCS`
        The best-fit WCS to the points given.
    """
    ...

def obsgeo_to_frame(obsgeo, obstime): # -> ITRS:
    """
    Convert a WCS obsgeo property into an `~.builtin_frames.ITRS` coordinate frame.

    Parameters
    ----------
    obsgeo : array-like
        A shape ``(6, )`` array representing ``OBSGEO-[XYZ], OBSGEO-[BLH]`` as
        returned by ``WCS.wcs.obsgeo``.

    obstime : time-like
        The time associated with the coordinate, will be passed to
        `~.builtin_frames.ITRS` as the obstime keyword.

    Returns
    -------
    `~.builtin_frames.ITRS`
        An `~.builtin_frames.ITRS` coordinate frame
        representing the coordinates.

    Notes
    -----

    The obsgeo array as accessed on a `.WCS` object is a length 6 numpy array
    where the first three elements are the coordinate in a cartesian
    representation and the second 3 are the coordinate in a spherical
    representation.

    This function priorities reading the cartesian coordinates, and will only
    read the spherical coordinates if the cartesian coordinates are either all
    zero or any of the cartesian coordinates are non-finite.

    In the case where both the spherical and cartesian coordinates have some
    non-finite values the spherical coordinates will be returned with the
    non-finite values included.

    """
    ...

