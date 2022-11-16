"""
This type stub file was generated by pyright.
"""

from astropy.utils.decorators import classproperty
from astropy.utils.state import ScienceState

"""
This module contains convenience functions for retrieving solar system
ephemerides from jplephem.
"""
__all__ = ["get_body", "get_moon", "get_body_barycentric", "get_body_barycentric_posvel", "solar_system_ephemeris"]
DEFAULT_JPL_EPHEMERIS = ...
BODY_NAME_TO_KERNEL_SPEC = ...
PLAN94_BODY_NAME_TO_PLANET_INDEX = ...
_EPHEMERIS_NOTE = ...
class solar_system_ephemeris(ScienceState):
    """Default ephemerides for calculating positions of Solar-System bodies.

    This can be one of the following:

    - 'builtin': polynomial approximations to the orbital elements.
    - 'dexxx[s]', for a JPL dynamical model, where xxx is the three digit
      version number (e.g. de430), and the 's' is optional to specify the
      'small' version of a kernel. The version number must correspond to an
      ephemeris file available at:
      https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/
    - 'jpl': Alias for the default JPL ephemeris (currently, 'de430').
    - URL: (str) The url to a SPK ephemeris in SPICE binary (.bsp) format.
    - PATH: (str) File path to a SPK ephemeris in SPICE binary (.bsp) format.
    - `None`: Ensure an Exception is raised without an explicit ephemeris.

    The default is 'builtin', which uses the ``epv00`` and ``plan94``
    routines from the ``erfa`` implementation of the Standards Of Fundamental
    Astronomy library.

    Notes
    -----
    Any file required will be downloaded (and cached) when the state is set.
    The default Satellite Planet Kernel (SPK) file from NASA JPL (de430) is
    ~120MB, and covers years ~1550-2650 CE [1]_.  The smaller de432s file is
    ~10MB, and covers years 1950-2050 [2]_ (and similarly for the newer de440
    and de440s).  Older versions of the JPL ephemerides (such as the widely
    used de200) can be used via their URL [3]_.

    .. [1] https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/aareadme_de430-de431.txt
    .. [2] https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/aareadme_de432s.txt
    .. [3] https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/
    """
    _value = ...
    _kernel = ...
    @classmethod
    def validate(cls, value): # -> str:
        ...
    
    @classmethod
    def get_kernel(cls, value): # -> None:
        ...
    
    @classproperty
    def kernel(cls): # -> None:
        ...
    
    @classproperty
    def bodies(cls): # -> tuple[str, ...] | None:
        ...
    


def get_body_barycentric_posvel(body, time, ephemeris=...): # -> tuple[CartesianRepresentation, CartesianRepresentation | Unbound] | CartesianRepresentation:
    """Calculate the barycentric position and velocity of a solar system body.

    Parameters
    ----------
    body : str or list of tuple
        The solar system body for which to calculate positions.  Can also be a
        kernel specifier (list of 2-tuples) if the ``ephemeris`` is a JPL
        kernel.
    time : `~astropy.time.Time`
        Time of observation.
    ephemeris : str, optional
        Ephemeris to use.  By default, use the one set with
        ``astropy.coordinates.solar_system_ephemeris.set``

    Returns
    -------
    position, velocity : tuple of `~astropy.coordinates.CartesianRepresentation`
        Tuple of barycentric (ICRS) position and velocity.

    See Also
    --------
    get_body_barycentric : to calculate position only.
        This is faster by about a factor two for JPL kernels, but has no
        speed advantage for the built-in ephemeris.

    Notes
    -----
    {_EPHEMERIS_NOTE}
    """
    ...

def get_body_barycentric(body, time, ephemeris=...): # -> tuple[CartesianRepresentation, CartesianRepresentation | Unbound] | CartesianRepresentation:
    """Calculate the barycentric position of a solar system body.

    Parameters
    ----------
    body : str or list of tuple
        The solar system body for which to calculate positions.  Can also be a
        kernel specifier (list of 2-tuples) if the ``ephemeris`` is a JPL
        kernel.
    time : `~astropy.time.Time`
        Time of observation.
    ephemeris : str, optional
        Ephemeris to use.  By default, use the one set with
        ``astropy.coordinates.solar_system_ephemeris.set``

    Returns
    -------
    position : `~astropy.coordinates.CartesianRepresentation`
        Barycentric (ICRS) position of the body in cartesian coordinates

    See Also
    --------
    get_body_barycentric_posvel : to calculate both position and velocity.

    Notes
    -----
    {_EPHEMERIS_NOTE}
    """
    ...

def get_body(body, time, location=..., ephemeris=...): # -> SkyCoord:
    """
    Get a `~astropy.coordinates.SkyCoord` for a solar system body as observed
    from a location on Earth in the `~astropy.coordinates.GCRS` reference
    system.

    Parameters
    ----------
    body : str or list of tuple
        The solar system body for which to calculate positions.  Can also be a
        kernel specifier (list of 2-tuples) if the ``ephemeris`` is a JPL
        kernel.
    time : `~astropy.time.Time`
        Time of observation.
    location : `~astropy.coordinates.EarthLocation`, optional
        Location of observer on the Earth.  If not given, will be taken from
        ``time`` (if not present, a geocentric observer will be assumed).
    ephemeris : str, optional
        Ephemeris to use.  If not given, use the one set with
        ``astropy.coordinates.solar_system_ephemeris.set`` (which is
        set to 'builtin' by default).

    Returns
    -------
    skycoord : `~astropy.coordinates.SkyCoord`
        GCRS Coordinate for the body

    Notes
    -----
    The coordinate returned is the apparent position, which is the position of
    the body at time *t* minus the light travel time from the *body* to the
    observing *location*.

    {_EPHEMERIS_NOTE}
    """
    ...

def get_moon(time, location=..., ephemeris=...): # -> SkyCoord:
    """
    Get a `~astropy.coordinates.SkyCoord` for the Earth's Moon as observed
    from a location on Earth in the `~astropy.coordinates.GCRS` reference
    system.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation
    location : `~astropy.coordinates.EarthLocation`
        Location of observer on the Earth. If none is supplied, taken from
        ``time`` (if not present, a geocentric observer will be assumed).
    ephemeris : str, optional
        Ephemeris to use.  If not given, use the one set with
        ``astropy.coordinates.solar_system_ephemeris.set`` (which is
        set to 'builtin' by default).

    Returns
    -------
    skycoord : `~astropy.coordinates.SkyCoord`
        GCRS Coordinate for the Moon

    Notes
    -----
    The coordinate returned is the apparent position, which is the position of
    the moon at time *t* minus the light travel time from the moon to the
    observing *location*.

    {_EPHEMERIS_NOTE}
    """
    ...

deprecation_msg = ...
