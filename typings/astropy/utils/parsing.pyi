"""
This type stub file was generated by pyright.
"""

import threading

"""
Wrappers for PLY to provide thread safety.
"""
__all__ = ['lex', 'ThreadSafeParser', 'yacc']
_TAB_HEADER = ...
_LOCK = threading.RLock()
def lex(lextab, package, reflags=...):
    """Create a lexer from local variables.

    It automatically compiles the lexer in optimized mode, writing to
    ``lextab`` in the same directory as the calling file.

    This function is thread-safe. The returned lexer is *not* thread-safe, but
    if it is used exclusively with a single parser returned by :func:`yacc`
    then it will be safe.

    It is only intended to work with lexers defined within the calling
    function, rather than at class or module scope.

    Parameters
    ----------
    lextab : str
        Name for the file to write with the generated tables, if it does not
        already exist (without ``.py`` suffix).
    package : str
        Name of a test package which should be run with pytest to regenerate
        the output file. This is inserted into a comment in the generated
        file.
    reflags : int
        Passed to ``ply.lex``.
    """
    ...

class ThreadSafeParser:
    """Wrap a parser produced by ``ply.yacc.yacc``.

    It provides a :meth:`parse` method that is thread-safe.
    """
    def __init__(self, parser) -> None:
        ...
    
    def parse(self, *args, **kwargs):
        """Run the wrapped parser, with a lock to ensure serialization."""
        ...
    


def yacc(tabmodule, package): # -> ThreadSafeParser:
    """Create a parser from local variables.

    It automatically compiles the parser in optimized mode, writing to
    ``tabmodule`` in the same directory as the calling file.

    This function is thread-safe, and the returned parser is also thread-safe,
    provided that it does not share a lexer with any other parser.

    It is only intended to work with parsers defined within the calling
    function, rather than at class or module scope.

    Parameters
    ----------
    tabmodule : str
        Name for the file to write with the generated tables, if it does not
        already exist (without ``.py`` suffix).
    package : str
        Name of a test package which should be run with pytest to regenerate
        the output file. This is inserted into a comment in the generated
        file.
    """
    ...
