"""
This type stub file was generated by pyright.
"""

__all__ = []
class HumanError(ValueError):
    ...


class CelestialError(ValueError):
    ...


def get_sign(dt): # -> Literal['capricorn', 'aquarius', 'pisces', 'aries', 'taurus', 'gemini', 'cancer', 'leo', 'virgo', 'libra', 'scorpio', 'sagittarius']:
    """
    """
    ...

_VALID_SIGNS = ...
_CONST_TO_SIGNS = ...
_ZODIAC = ...
def horoscope(birthday, corrected=..., chinese=...): # -> None:
    """
    Enter your birthday as an `astropy.time.Time` object and
    receive a mystical horoscope about things to come.

    Parameters
    ----------
    birthday : `astropy.time.Time` or str
        Your birthday as a `datetime.datetime` or `astropy.time.Time` object
        or "YYYY-MM-DD"string.
    corrected : bool
        Whether to account for the precession of the Earth instead of using the
        ancient Greek dates for the signs.  After all, you do want your *real*
        horoscope, not a cheap inaccurate approximation, right?

    chinese : bool
        Chinese annual zodiac wisdom instead of Western one.

    Returns
    -------
    Infinite wisdom, condensed into astrologically precise prose.

    Notes
    -----
    This function was implemented on April 1.  Take note of that date.
    """
    ...

def inject_horoscope(): # -> None:
    ...

