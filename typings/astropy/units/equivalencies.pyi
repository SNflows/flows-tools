"""
This type stub file was generated by pyright.
"""

from collections import UserList

"""A set of standard astronomical equivalencies."""
__all__ = ['parallax', 'spectral', 'spectral_density', 'doppler_radio', 'doppler_optical', 'doppler_relativistic', 'doppler_redshift', 'mass_energy', 'brightness_temperature', 'thermodynamic_temperature', 'beam_angular_area', 'dimensionless_angles', 'logarithmic', 'temperature', 'temperature_energy', 'molar_mass_amu', 'pixel_scale', 'plate_scale', "Equivalency"]
class Equivalency(UserList):
    """
    A container for a units equivalency.

    Attributes
    ----------
    name: `str`
        The name of the equivalency.
    kwargs: `dict`
        Any positional or keyword arguments used to make the equivalency.
    """
    def __init__(self, equiv_list, name=..., kwargs=...) -> None:
        ...
    
    def __add__(self, other): # -> Equivalency:
        ...
    
    def __eq__(self, other) -> bool:
        ...
    


def dimensionless_angles(): # -> Equivalency:
    """Allow angles to be equivalent to dimensionless (with 1 rad = 1 m/m = 1).

    It is special compared to other equivalency pairs in that it
    allows this independent of the power to which the angle is raised,
    and independent of whether it is part of a more complicated unit.
    """
    ...

def logarithmic(): # -> Equivalency:
    """Allow logarithmic units to be converted to dimensionless fractions"""
    ...

def parallax(): # -> Equivalency:
    """
    Returns a list of equivalence pairs that handle the conversion
    between parallax angle and distance.
    """
    ...

def spectral(): # -> Equivalency:
    """
    Returns a list of equivalence pairs that handle spectral
    wavelength, wave number, frequency, and energy equivalencies.

    Allows conversions between wavelength units, wave number units,
    frequency units, and energy units as they relate to light.

    There are two types of wave number:

        * spectroscopic - :math:`1 / \\lambda` (per meter)
        * angular - :math:`2 \\pi / \\lambda` (radian per meter)

    """
    ...

def spectral_density(wav, factor=...): # -> Equivalency:
    """
    Returns a list of equivalence pairs that handle spectral density
    with regard to wavelength and frequency.

    Parameters
    ----------
    wav : `~astropy.units.Quantity`
        `~astropy.units.Quantity` associated with values being converted
        (e.g., wavelength or frequency).

    Notes
    -----
    The ``factor`` argument is left for backward-compatibility with the syntax
    ``spectral_density(unit, factor)`` but users are encouraged to use
    ``spectral_density(factor * unit)`` instead.

    """
    ...

def doppler_radio(rest): # -> Equivalency:
    r"""
    Return the equivalency pairs for the radio convention for velocity.

    The radio convention for the relation between velocity and frequency is:

    :math:`V = c \frac{f_0 - f}{f_0}  ;  f(V) = f_0 ( 1 - V/c )`

    Parameters
    ----------
    rest : `~astropy.units.Quantity`
        Any quantity supported by the standard spectral equivalencies
        (wavelength, energy, frequency, wave number).

    References
    ----------
    `NRAO site defining the conventions <https://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html>`_

    Examples
    --------
    >>> import astropy.units as u
    >>> CO_restfreq = 115.27120*u.GHz  # rest frequency of 12 CO 1-0 in GHz
    >>> radio_CO_equiv = u.doppler_radio(CO_restfreq)
    >>> measured_freq = 115.2832*u.GHz
    >>> radio_velocity = measured_freq.to(u.km/u.s, equivalencies=radio_CO_equiv)
    >>> radio_velocity  # doctest: +FLOAT_CMP
    <Quantity -31.209092088877583 km / s>
    """
    ...

def doppler_optical(rest): # -> Equivalency:
    r"""
    Return the equivalency pairs for the optical convention for velocity.

    The optical convention for the relation between velocity and frequency is:

    :math:`V = c \frac{f_0 - f}{f  }  ;  f(V) = f_0 ( 1 + V/c )^{-1}`

    Parameters
    ----------
    rest : `~astropy.units.Quantity`
        Any quantity supported by the standard spectral equivalencies
        (wavelength, energy, frequency, wave number).

    References
    ----------
    `NRAO site defining the conventions <https://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html>`_

    Examples
    --------
    >>> import astropy.units as u
    >>> CO_restfreq = 115.27120*u.GHz  # rest frequency of 12 CO 1-0 in GHz
    >>> optical_CO_equiv = u.doppler_optical(CO_restfreq)
    >>> measured_freq = 115.2832*u.GHz
    >>> optical_velocity = measured_freq.to(u.km/u.s, equivalencies=optical_CO_equiv)
    >>> optical_velocity  # doctest: +FLOAT_CMP
    <Quantity -31.20584348799674 km / s>
    """
    ...

def doppler_relativistic(rest): # -> Equivalency:
    r"""
    Return the equivalency pairs for the relativistic convention for velocity.

    The full relativistic convention for the relation between velocity and frequency is:

    :math:`V = c \frac{f_0^2 - f^2}{f_0^2 + f^2} ;  f(V) = f_0 \frac{\left(1 - (V/c)^2\right)^{1/2}}{(1+V/c)}`

    Parameters
    ----------
    rest : `~astropy.units.Quantity`
        Any quantity supported by the standard spectral equivalencies
        (wavelength, energy, frequency, wave number).

    References
    ----------
    `NRAO site defining the conventions <https://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html>`_

    Examples
    --------
    >>> import astropy.units as u
    >>> CO_restfreq = 115.27120*u.GHz  # rest frequency of 12 CO 1-0 in GHz
    >>> relativistic_CO_equiv = u.doppler_relativistic(CO_restfreq)
    >>> measured_freq = 115.2832*u.GHz
    >>> relativistic_velocity = measured_freq.to(u.km/u.s, equivalencies=relativistic_CO_equiv)
    >>> relativistic_velocity  # doctest: +FLOAT_CMP
    <Quantity -31.207467619351537 km / s>
    >>> measured_velocity = 1250 * u.km/u.s
    >>> relativistic_frequency = measured_velocity.to(u.GHz, equivalencies=relativistic_CO_equiv)
    >>> relativistic_frequency  # doctest: +FLOAT_CMP
    <Quantity 114.79156866993588 GHz>
    >>> relativistic_wavelength = measured_velocity.to(u.mm, equivalencies=relativistic_CO_equiv)
    >>> relativistic_wavelength  # doctest: +FLOAT_CMP
    <Quantity 2.6116243681798923 mm>
    """
    ...

def doppler_redshift(): # -> Equivalency:
    """
    Returns the equivalence between Doppler redshift (unitless) and radial velocity.

    .. note::

        This equivalency is not compatible with cosmological
        redshift in `astropy.cosmology.units`.

    """
    ...

def molar_mass_amu(): # -> Equivalency:
    """
    Returns the equivalence between amu and molar mass.
    """
    ...

def mass_energy(): # -> Equivalency:
    """
    Returns a list of equivalence pairs that handle the conversion
    between mass and energy.
    """
    ...

def brightness_temperature(frequency, beam_area=...): # -> Equivalency:
    r"""
    Defines the conversion between Jy/sr and "brightness temperature",
    :math:`T_B`, in Kelvins.  The brightness temperature is a unit very
    commonly used in radio astronomy.  See, e.g., "Tools of Radio Astronomy"
    (Wilson 2009) eqn 8.16 and eqn 8.19 (these pages are available on `google
    books
    <https://books.google.com/books?id=9KHw6R8rQEMC&pg=PA179&source=gbs_toc_r&cad=4#v=onepage&q&f=false>`__).

    :math:`T_B \equiv S_\nu / \left(2 k \nu^2 / c^2 \right)`

    If the input is in Jy/beam or Jy (assuming it came from a single beam), the
    beam area is essential for this computation: the brightness temperature is
    inversely proportional to the beam area.

    Parameters
    ----------
    frequency : `~astropy.units.Quantity`
        The observed ``spectral`` equivalent `~astropy.units.Unit` (e.g.,
        frequency or wavelength).  The variable is named 'frequency' because it
        is more commonly used in radio astronomy.
        BACKWARD COMPATIBILITY NOTE: previous versions of the brightness
        temperature equivalency used the keyword ``disp``, which is no longer
        supported.
    beam_area : `~astropy.units.Quantity` ['solid angle']
        Beam area in angular units, i.e. steradian equivalent

    Examples
    --------
    Arecibo C-band beam::

        >>> import numpy as np
        >>> from astropy import units as u
        >>> beam_sigma = 50*u.arcsec
        >>> beam_area = 2*np.pi*(beam_sigma)**2
        >>> freq = 5*u.GHz
        >>> equiv = u.brightness_temperature(freq)
        >>> (1*u.Jy/beam_area).to(u.K, equivalencies=equiv)  # doctest: +FLOAT_CMP
        <Quantity 3.526295144567176 K>

    VLA synthetic beam::

        >>> bmaj = 15*u.arcsec
        >>> bmin = 15*u.arcsec
        >>> fwhm_to_sigma = 1./(8*np.log(2))**0.5
        >>> beam_area = 2.*np.pi*(bmaj*bmin*fwhm_to_sigma**2)
        >>> freq = 5*u.GHz
        >>> equiv = u.brightness_temperature(freq)
        >>> (u.Jy/beam_area).to(u.K, equivalencies=equiv)  # doctest: +FLOAT_CMP
        <Quantity 217.2658703625732 K>

    Any generic surface brightness:

        >>> surf_brightness = 1e6*u.MJy/u.sr
        >>> surf_brightness.to(u.K, equivalencies=u.brightness_temperature(500*u.GHz)) # doctest: +FLOAT_CMP
        <Quantity 130.1931904778803 K>
    """
    ...

def beam_angular_area(beam_area): # -> Equivalency:
    """
    Convert between the ``beam`` unit, which is commonly used to express the area
    of a radio telescope resolution element, and an area on the sky.
    This equivalency also supports direct conversion between ``Jy/beam`` and
    ``Jy/steradian`` units, since that is a common operation.

    Parameters
    ----------
    beam_area : unit-like
        The area of the beam in angular area units (e.g., steradians)
        Must have angular area equivalent units.
    """
    ...

def thermodynamic_temperature(frequency, T_cmb=...): # -> Equivalency:
    r"""Defines the conversion between Jy/sr and "thermodynamic temperature",
    :math:`T_{CMB}`, in Kelvins.  The thermodynamic temperature is a unit very
    commonly used in cosmology. See eqn 8 in [1]

    :math:`K_{CMB} \equiv I_\nu / \left(2 k \nu^2 / c^2  f(\nu) \right)`

    with :math:`f(\nu) = \frac{ x^2 e^x}{(e^x - 1 )^2}`
    where :math:`x = h \nu / k T`

    Parameters
    ----------
    frequency : `~astropy.units.Quantity`
        The observed `spectral` equivalent `~astropy.units.Unit` (e.g.,
        frequency or wavelength). Must have spectral units.
    T_cmb :  `~astropy.units.Quantity` ['temperature'] or None
        The CMB temperature at z=0.  If `None`, the default cosmology will be
        used to get this temperature. Must have units of temperature.

    Notes
    -----
    For broad band receivers, this conversion do not hold
    as it highly depends on the frequency

    References
    ----------
    .. [1] Planck 2013 results. IX. HFI spectral response
       https://arxiv.org/abs/1303.5070

    Examples
    --------
    Planck HFI 143 GHz::

        >>> from astropy import units as u
        >>> from astropy.cosmology import Planck15
        >>> freq = 143 * u.GHz
        >>> equiv = u.thermodynamic_temperature(freq, Planck15.Tcmb0)
        >>> (1. * u.mK).to(u.MJy / u.sr, equivalencies=equiv)  # doctest: +FLOAT_CMP
        <Quantity 0.37993172 MJy / sr>

    """
    ...

def temperature(): # -> Equivalency:
    """Convert between Kelvin, Celsius, Rankine and Fahrenheit here because
    Unit and CompositeUnit cannot do addition or subtraction properly.
    """
    ...

def temperature_energy(): # -> Equivalency:
    """Convert between Kelvin and keV(eV) to an equivalent amount."""
    ...

def assert_is_spectral_unit(value): # -> None:
    ...

def pixel_scale(pixscale): # -> Equivalency:
    """
    Convert between pixel distances (in units of ``pix``) and other units,
    given a particular ``pixscale``.

    Parameters
    ----------
    pixscale : `~astropy.units.Quantity`
        The pixel scale either in units of <unit>/pixel or pixel/<unit>.
    """
    ...

def plate_scale(platescale): # -> Equivalency:
    """
    Convert between lengths (to be interpreted as lengths in the focal plane)
    and angular units with a specified ``platescale``.

    Parameters
    ----------
    platescale : `~astropy.units.Quantity`
        The pixel scale either in units of distance/pixel or distance/angle.
    """
    ...

def __getattr__(attr): # -> (H0: Unknown | None = None) -> Equivalency:
    ...
