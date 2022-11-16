"""
This type stub file was generated by pyright.
"""

from functools import lru_cache
from typing import Optional, Union
from enum import Enum

"""
Get FLOWS specific catalogs.
"""

class OutputFormat(Enum):
    table = ...
    dictionary = ...
    json = ...

def create_catalog_table(jsn: dict) -> dict: ...
@lru_cache(maxsize=10)
def get_catalog(target: Union[int, str], radius: Optional[float] = ..., output: str = ...) -> dict:
    """

    Parameters:
        target (int or str):
        radius (float, optional): Radius around target in degrees to return targets for.
        output (str, optional): Desired output format. Choices are 'table', 'dict', 'json'.
            Default='table'.

    Returns:
        dict: Dictionary with three members:
            - 'target': Information about target.
            - 'references': Table with information about reference stars close to target.
            - 'avoid': Table with stars close to target which should be avoided in FOV selection.
    """
    ...

def get_catalog_missing():  # -> Any:
    """
    Get missing catalogs
    Returns: json of missing catalogs
    """
    ...