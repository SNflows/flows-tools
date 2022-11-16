"""
This type stub file was generated by pyright.
"""

"""
A collection of functions for checking various XML-related strings for
standards compliance.
"""
def check_id(ID): # -> bool:
    """
    Returns `True` if *ID* is a valid XML ID.
    """
    ...

def fix_id(ID): # -> str:
    """
    Given an arbitrary string, create one that can be used as an xml
    id.  This is rather simplistic at the moment, since it just
    replaces non-valid characters with underscores.
    """
    ...

_token_regex = ...
def check_token(token): # -> bool:
    """
    Returns `True` if *token* is a valid XML token, as defined by XML
    Schema Part 2.
    """
    ...

def check_mime_content_type(content_type): # -> bool:
    """
    Returns `True` if *content_type* is a valid MIME content type
    (syntactically at least), as defined by RFC 2045.
    """
    ...

def check_anyuri(uri): # -> bool:
    """
    Returns `True` if *uri* is a valid URI as defined in RFC 2396.
    """
    ...

