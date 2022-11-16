from astropy.utils.exceptions import AstropyUserWarning

class CRTFRegionParserWarning(AstropyUserWarning): ...
class CRTFRegionParserError(ValueError): ...

# Valid symbols in CRTF
valid_symbols: dict[str, str]