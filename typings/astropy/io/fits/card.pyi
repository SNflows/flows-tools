"""
This type stub file was generated by pyright.
"""

from .verify import _Verify

__all__ = ['Card', 'Undefined']
FIX_FP_TABLE = ...
FIX_FP_TABLE2 = ...
CARD_LENGTH = ...
BLANK_CARD = ...
KEYWORD_LENGTH = ...
VALUE_INDICATOR = ...
VALUE_INDICATOR_LEN = ...
HIERARCH_VALUE_INDICATOR = ...
class Undefined:
    """Undefined value."""
    def __init__(self) -> None:
        ...
    


UNDEFINED = ...
class Card(_Verify):
    length = ...
    _keywd_FSC_RE = ...
    _keywd_hierarch_RE = ...
    _digits_FSC = ...
    _digits_NFSC = ...
    _numr_FSC = ...
    _numr_NFSC = ...
    _number_FSC_RE = ...
    _number_NFSC_RE = ...
    _strg = ...
    _comm_field = ...
    _strg_comment_RE = ...
    _ascii_text_re = ...
    _value_FSC_RE = ...
    _value_NFSC_RE = ...
    _rvkc_identifier = ...
    _rvkc_field = ...
    _rvkc_field_specifier_s = ...
    _rvkc_field_specifier_val = ...
    _rvkc_keyword_val = ...
    _rvkc_keyword_val_comm = ...
    _rvkc_field_specifier_val_RE = ...
    _rvkc_keyword_name_RE = ...
    _rvkc_keyword_val_comm_RE = ...
    _commentary_keywords = ...
    _special_keywords = ...
    _value_indicator = ...
    def __init__(self, keyword=..., value=..., comment=..., **kwargs) -> None:
        ...
    
    def __repr__(self): # -> str:
        ...
    
    def __str__(self) -> str:
        ...
    
    def __len__(self): # -> Literal[3]:
        ...
    
    def __getitem__(self, index): # -> str | float | Any | Undefined | bool | float32 | int | complex | LiteralString | None:
        ...
    
    @property
    def keyword(self): # -> str | None:
        """Returns the keyword name parsed from the card image."""
        ...
    
    @keyword.setter
    def keyword(self, keyword): # -> None:
        """Set the key attribute; once set it cannot be modified."""
        ...
    
    @property
    def value(self): # -> float | str | Any | Undefined | bool | float32 | int | complex | None:
        """The value associated with the keyword stored in this card."""
        ...
    
    @value.setter
    def value(self, value): # -> None:
        ...
    
    @value.deleter
    def value(self): # -> None:
        ...
    
    @property
    def rawkeyword(self): # -> str | None:
        """On record-valued keyword cards this is the name of the standard <= 8
        character FITS keyword that this RVKC is stored in.  Otherwise it is
        the card's normal keyword.
        """
        ...
    
    @property
    def rawvalue(self): # -> str | float | Any | Undefined | bool | float32 | int | complex | None:
        """On record-valued keyword cards this is the raw string value in
        the ``<field-specifier>: <value>`` format stored in the card in order
        to represent a RVKC.  Otherwise it is the card's normal value.
        """
        ...
    
    @property
    def comment(self): # -> str | Any | LiteralString:
        """Get the comment attribute from the card image if not already set."""
        ...
    
    @comment.setter
    def comment(self, comment): # -> None:
        ...
    
    @comment.deleter
    def comment(self): # -> None:
        ...
    
    @property
    def field_specifier(self): # -> None:
        """
        The field-specifier of record-valued keyword cards; always `None` on
        normal cards.
        """
        ...
    
    @field_specifier.setter
    def field_specifier(self, field_specifier): # -> None:
        ...
    
    @field_specifier.deleter
    def field_specifier(self):
        ...
    
    @property
    def image(self):
        """
        The card "image", that is, the 80 byte character string that represents
        this card in an actual FITS header.
        """
        ...
    
    @property
    def is_blank(self): # -> bool:
        """
        `True` if the card is completely blank--that is, it has no keyword,
        value, or comment.  It appears in the header as 80 spaces.

        Returns `False` otherwise.
        """
        ...
    
    @classmethod
    def fromstring(cls, image): # -> Self@Card:
        """
        Construct a `Card` object from a (raw) string. It will pad the string
        if it is not the length of a card image (80 columns).  If the card
        image is longer than 80 columns, assume it contains ``CONTINUE``
        card(s).
        """
        ...
    
    @classmethod
    def normalize_keyword(cls, keyword): # -> str:
        """
        `classmethod` to convert a keyword value that may contain a
        field-specifier to uppercase.  The effect is to raise the key to
        uppercase and leave the field specifier in its original case.

        Parameters
        ----------
        keyword : or str
            A keyword value or a ``keyword.field-specifier`` value
        """
        ...
    


