"""
This type stub file was generated by pyright.
"""

DMS_RE = ...
HMS_RE = ...
DDEC_RE = ...
DMIN_RE = ...
DSEC_RE = ...
SCAL_RE = ...
CUSTOM_UNITS = ...
class BaseFormatterLocator:
    """
    A joint formatter/locator
    """
    def __init__(self, values=..., number=..., spacing=..., format=..., unit=..., format_unit=...) -> None:
        ...
    
    @property
    def values(self): # -> Quantity | None:
        ...
    
    @values.setter
    def values(self, values): # -> None:
        ...
    
    @property
    def number(self): # -> None:
        ...
    
    @number.setter
    def number(self, number): # -> None:
        ...
    
    @property
    def spacing(self): # -> None:
        ...
    
    @spacing.setter
    def spacing(self, spacing): # -> None:
        ...
    
    def minor_locator(self, spacing, frequency, value_min, value_max):
        ...
    
    @property
    def format_unit(self): # -> Unit | None:
        ...
    
    @format_unit.setter
    def format_unit(self, unit): # -> None:
        ...
    


class AngleFormatterLocator(BaseFormatterLocator):
    """
    A joint formatter/locator
    """
    def __init__(self, values=..., number=..., spacing=..., format=..., unit=..., decimal=..., format_unit=..., show_decimal_unit=...) -> None:
        ...
    
    @property
    def decimal(self): # -> bool | None:
        ...
    
    @decimal.setter
    def decimal(self, value): # -> None:
        ...
    
    @property
    def spacing(self): # -> Quantity:
        ...
    
    @spacing.setter
    def spacing(self, spacing): # -> None:
        ...
    
    @property
    def sep(self): # -> None:
        ...
    
    @sep.setter
    def sep(self, separator): # -> None:
        ...
    
    @property
    def format(self):
        ...
    
    @format.setter
    def format(self, value): # -> None:
        ...
    
    @property
    def base_spacing(self): # -> int:
        ...
    
    def locator(self, value_min, value_max): # -> tuple[Quantity, Unknown] | tuple[Unknown, Unknown]:
        ...
    
    def formatter(self, values, spacing, format=...): # -> Any | list[Unknown]:
        ...
    


class ScalarFormatterLocator(BaseFormatterLocator):
    """
    A joint formatter/locator
    """
    def __init__(self, values=..., number=..., spacing=..., format=..., unit=..., format_unit=...) -> None:
        ...
    
    @property
    def spacing(self): # -> Quantity:
        ...
    
    @spacing.setter
    def spacing(self, spacing): # -> None:
        ...
    
    @property
    def format(self):
        ...
    
    @format.setter
    def format(self, value): # -> None:
        ...
    
    @property
    def base_spacing(self):
        ...
    
    def locator(self, value_min, value_max): # -> tuple[Quantity, Unknown] | tuple[Unknown, Unknown]:
        ...
    
    def formatter(self, values, spacing, format=...): # -> list[Unknown] | list[str]:
        ...
    

