"""
This type stub file was generated by pyright.
"""

from . import base

"""
Handles the "LaTeX" unit format.
"""
class Latex(base.Base):
    """
    Output LaTeX to display the unit based on IAU style guidelines.

    Attempts to follow the `IAU Style Manual
    <https://www.iau.org/static/publications/stylemanual1989.pdf>`_.
    """
    @classmethod
    def to_string(cls, unit): # -> str:
        ...
    
    @classmethod
    def format_exponential_notation(cls, val, format_spec=...): # -> LiteralString | Literal['{\\rm NaN}', '\\infty', '-\\infty']:
        """
        Formats a value in exponential notation for LaTeX.

        Parameters
        ----------
        val : number
            The value to be formatted

        format_spec : str, optional
            Format used to split up mantissa and exponent

        Returns
        -------
        latex_string : str
            The value in exponential notation in a format suitable for LaTeX.
        """
        ...
    


class LatexInline(Latex):
    """
    Output LaTeX to display the unit based on IAU style guidelines with negative
    powers.

    Attempts to follow the `IAU Style Manual
    <https://www.iau.org/static/publications/stylemanual1989.pdf>`_ and the
    `ApJ and AJ style guide
    <https://journals.aas.org/manuscript-preparation/>`_.
    """
    name = ...


