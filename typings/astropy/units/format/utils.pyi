"""
This type stub file was generated by pyright.
"""

"""
Utilities shared by the different formats.
"""
def get_grouped_by_powers(bases, powers): # -> tuple[list[Unknown], list[Unknown]]:
    """
    Groups the powers and bases in the given
    `~astropy.units.CompositeUnit` into positive powers and
    negative powers for easy display on either side of a solidus.

    Parameters
    ----------
    bases : list of `astropy.units.UnitBase` instances

    powers : list of int

    Returns
    -------
    positives, negatives : tuple of lists
       Each element in each list is tuple of the form (*base*,
       *power*).  The negatives have the sign of their power reversed
       (i.e. the powers are all positive).
    """
    ...

def split_mantissa_exponent(v, format_spec=...): # -> tuple[str, str]:
    """
    Given a number, split it into its mantissa and base 10 exponent
    parts, each as strings.  If the exponent is too small, it may be
    returned as the empty string.

    Parameters
    ----------
    v : float

    format_spec : str, optional
        Number representation formatting string

    Returns
    -------
    mantissa, exponent : tuple of strings
    """
    ...

def decompose_to_known_units(unit, func): # -> Unit | NamedUnit:
    """
    Partially decomposes a unit so it is only composed of units that
    are "known" to a given format.

    Parameters
    ----------
    unit : `~astropy.units.UnitBase` instance

    func : callable
        This function will be called to determine if a given unit is
        "known".  If the unit is not known, this function should raise a
        `ValueError`.

    Returns
    -------
    unit : `~astropy.units.UnitBase` instance
        A flattened unit.
    """
    ...

def format_power(power): # -> str:
    """
    Converts a value for a power (which may be floating point or a
    `fractions.Fraction` object), into a string looking like either
    an integer or a fraction, if the power is close to that.
    """
    ...

def did_you_mean_units(s, all_units, deprecated_units, format_decomposed): # -> str:
    """
    A wrapper around `astropy.utils.misc.did_you_mean` that deals with
    the display of deprecated units.

    Parameters
    ----------
    s : str
        The invalid unit string

    all_units : dict
        A mapping from valid unit names to unit objects.

    deprecated_units : sequence
        The deprecated unit names

    format_decomposed : callable
        A function to turn a decomposed version of the unit into a
        string.  Should return `None` if not possible

    Returns
    -------
    msg : str
        A string message with a list of alternatives, or the empty
        string.
    """
    ...

def unit_deprecation_warning(s, unit, standard_name, format_decomposed): # -> None:
    """
    Raises a UnitsWarning about a deprecated unit in a given format.
    Suggests a decomposed alternative if one is available.

    Parameters
    ----------
    s : str
        The deprecated unit name.

    unit : astropy.units.core.UnitBase
        The unit object.

    standard_name : str
        The name of the format for which the unit is deprecated.

    format_decomposed : callable
        A function to turn a decomposed version of the unit into a
        string.  Should return `None` if not possible
    """
    ...

