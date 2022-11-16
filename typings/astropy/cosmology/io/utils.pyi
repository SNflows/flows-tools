"""
This type stub file was generated by pyright.
"""

FULLQUALNAME_SUBSTITUTIONS = ...
def convert_parameter_to_column(parameter, value, meta=...): # -> Column:
    """Convert a |Cosmology| Parameter to a Table |Column|.

    Parameters
    ----------
    parameter : `astropy.cosmology.parameter.Parameter`
    value : Any
    meta : dict or None, optional
        Information from the Cosmology's metadata.

    Returns
    -------
    `astropy.table.Column`
    """
    ...

def convert_parameter_to_model_parameter(parameter, value, meta=...): # -> Parameter:
    """Convert a Cosmology Parameter to a Model Parameter.

    Parameters
    ----------
    parameter : `astropy.cosmology.parameter.Parameter`
    value : Any
    meta : dict or None, optional
        Information from the Cosmology's metadata.
        This function will use any of: 'getter', 'setter', 'fixed', 'tied',
        'min', 'max', 'bounds', 'prior', 'posterior'.

    Returns
    -------
    `astropy.modeling.Parameter`
    """
    ...
