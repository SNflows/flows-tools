"""
This type stub file was generated by pyright.
"""

compiler = ...
BOMS = ...
BOM_LIST = ...
BOM_SET = ...
def match_utf8(encoding): # -> bool:
    ...

squot = ...
dquot = ...
noquot = ...
wspace_plus = ...
tsquot = ...
tdquot = ...
MISSING = ...
__all__ = ('DEFAULT_INDENT_TYPE', 'DEFAULT_INTERPOLATION', 'ConfigObjError', 'NestingError', 'ParseError', 'DuplicateError', 'ConfigspecError', 'ConfigObj', 'SimpleVal', 'InterpolationError', 'InterpolationLoopError', 'MissingInterpolationOption', 'RepeatSectionError', 'ReloadError', 'UnreprError', 'UnknownType', 'flatten_errors', 'get_extra_values')
DEFAULT_INTERPOLATION = ...
DEFAULT_INDENT_TYPE = ...
MAX_INTERPOL_DEPTH = ...
OPTION_DEFAULTS = ...
def getObj(s):
    ...

class UnknownType(Exception):
    ...


class Builder:
    def build(self, o):
        ...
    
    def build_List(self, o): # -> list[Unknown]:
        ...
    
    def build_Const(self, o):
        ...
    
    def build_Dict(self, o): # -> dict[Unknown, Unknown]:
        ...
    
    def build_Tuple(self, o): # -> tuple[Unknown, ...]:
        ...
    
    def build_Name(self, o): # -> bool | None:
        ...
    
    def build_Add(self, o): # -> complex:
        ...
    
    def build_Getattr(self, o): # -> Any:
        ...
    
    def build_UnarySub(self, o):
        ...
    
    def build_UnaryAdd(self, o):
        ...
    


_builder = ...
def unrepr(s): # -> Any:
    ...

class ConfigObjError(SyntaxError):
    """
    This is the base class for all errors that ConfigObj raises.
    It is a subclass of SyntaxError.
    """
    def __init__(self, message=..., line_number=..., line=...) -> None:
        ...
    


class NestingError(ConfigObjError):
    """
    This error indicates a level of nesting that doesn't match.
    """
    ...


class ParseError(ConfigObjError):
    """
    This error indicates that a line is badly written.
    It is neither a valid ``key = value`` line,
    nor a valid section marker line.
    """
    ...


class ReloadError(IOError):
    """
    A 'reload' operation failed.
    This exception is a subclass of ``IOError``.
    """
    def __init__(self) -> None:
        ...
    


class DuplicateError(ConfigObjError):
    """
    The keyword or section specified already exists.
    """
    ...


class ConfigspecError(ConfigObjError):
    """
    An error occured whilst parsing a configspec.
    """
    ...


class InterpolationError(ConfigObjError):
    """Base class for the two interpolation errors."""
    ...


class InterpolationLoopError(InterpolationError):
    """Maximum interpolation depth exceeded in string interpolation."""
    def __init__(self, option) -> None:
        ...
    


class RepeatSectionError(ConfigObjError):
    """
    This error indicates additional sections in a section with a
    ``__many__`` (repeated) section.
    """
    ...


class MissingInterpolationOption(InterpolationError):
    """A value specified for interpolation was missing."""
    def __init__(self, option) -> None:
        ...
    


class UnreprError(ConfigObjError):
    """An error parsing in unrepr mode."""
    ...


class InterpolationEngine:
    """
    A helper class to help perform string interpolation.

    This class is an abstract base class; its descendants perform
    the actual work.
    """
    _KEYCRE = ...
    _cookie = ...
    def __init__(self, section) -> None:
        ...
    
    def interpolate(self, key, value): # -> LiteralString:
        ...
    


class ConfigParserInterpolation(InterpolationEngine):
    """Behaves like ConfigParser."""
    _cookie = ...
    _KEYCRE = ...


class TemplateInterpolation(InterpolationEngine):
    """Behaves like string.Template."""
    _cookie = ...
    _delimiter = ...
    _KEYCRE = ...


interpolation_engines = ...
def __newobj__(cls, *args):
    ...

class Section(dict):
    """
    A dictionary-like object that represents a section in a config file.

    It does string interpolation if the 'interpolation' attribute
    of the 'main' object is set to True.

    Interpolation is tried first from this object, then from the 'DEFAULT'
    section of this object, next from the parent and its 'DEFAULT' section,
    and so on until the main object is reached.

    A Section will behave like an ordered dictionary - following the
    order of the ``scalars`` and ``sections`` attributes.
    You can use this to change the order of members.

    Iteration follows the order: scalars, then sections.
    """
    def __setstate__(self, state): # -> None:
        ...
    
    def __reduce__(self): # -> tuple[(cls: Unknown, *args: Unknown) -> Unknown, tuple[Type[Self@Section]], tuple[dict[Unknown, str | Unknown | list[str | Unknown] | list[Unknown]], dict[str, Any]]]:
        ...
    
    def __init__(self, parent, depth, main, indict=..., name=...) -> None:
        """
        * parent is the section above
        * depth is the depth level of this section
        * main is the main ConfigObj
        * indict is a dictionary to initialise the section with
        """
        ...
    
    def __getitem__(self, key): # -> str | LiteralString | list[str | Unknown] | list[Unknown]:
        """Fetch the item and do string interpolation."""
        ...
    
    def __setitem__(self, key, value, unrepr=...): # -> None:
        """
        Correctly set a value.

        Making dictionary values Section instances.
        (We have to special case 'Section' instances - which are also dicts)

        Keys must be strings.
        Values need only be strings (or lists of strings) if
        ``main.stringify`` is set.

        ``unrepr`` must be set when setting a value to a dictionary, without
        creating a new sub-section.
        """
        ...
    
    def __delitem__(self, key): # -> None:
        """Remove items from the sequence when deleting."""
        ...
    
    def get(self, key, default=...): # -> str | LiteralString | list[str | Unknown] | list[Unknown] | None:
        """A version of ``get`` that doesn't bypass string interpolation."""
        ...
    
    def update(self, indict): # -> None:
        """
        A version of update that uses our ``__setitem__``.
        """
        ...
    
    def pop(self, key, default=...): # -> str | LiteralString | list[str | Unknown] | list[Unknown] | object:
        """
        'D.pop(k[,d]) -> v, remove specified key and return the corresponding value.
        If key is not found, d is returned if given, otherwise KeyError is raised'
        """
        ...
    
    def popitem(self): # -> tuple[Unknown, str | Unknown | LiteralString | list[str | Unknown] | list[Unknown]]:
        """Pops the first (key,val)"""
        ...
    
    def clear(self): # -> None:
        """
        A version of clear that also affects scalars/sections
        Also clears comments and configspec.

        Leaves other attributes alone :
            depth/main/parent are not affected
        """
        ...
    
    def setdefault(self, key, default=...): # -> str | LiteralString | list[str | Unknown] | list[Unknown]:
        """A version of setdefault that sets sequence if appropriate."""
        ...
    
    def items(self): # -> list[tuple[Unknown, str | Unknown | list[str | Unknown] | list[Unknown]]]:
        """D.items() -> list of D's (key, value) pairs, as 2-tuples"""
        ...
    
    def keys(self): # -> list[Unknown]:
        """D.keys() -> list of D's keys"""
        ...
    
    def values(self): # -> list[str | Unknown | list[str | Unknown] | list[Unknown]]:
        """D.values() -> list of D's values"""
        ...
    
    def iteritems(self): # -> Iterator[tuple[Unknown, str | Unknown | list[str | Unknown] | list[Unknown]]]:
        """D.iteritems() -> an iterator over the (key, value) items of D"""
        ...
    
    def iterkeys(self): # -> Iterator[Unknown]:
        """D.iterkeys() -> an iterator over the keys of D"""
        ...
    
    __iter__ = ...
    def itervalues(self): # -> Iterator[str | Unknown | list[str | Unknown] | list[Unknown]]:
        """D.itervalues() -> an iterator over the values of D"""
        ...
    
    def __repr__(self): # -> str:
        """x.__repr__() <==> repr(x)"""
        ...
    
    __str__ = ...
    def dict(self): # -> dict[Unknown, Unknown]:
        """
        Return a deepcopy of self as a dictionary.

        All members that are ``Section`` instances are recursively turned to
        ordinary dictionaries - by calling their ``dict`` method.

        >>> n = a.dict()
        >>> n == a
        1
        >>> n is a
        0
        """
        ...
    
    def merge(self, indict): # -> None:
        """
        A recursive update - useful for merging config files.

        >>> a = '''[section1]
        ...     option1 = True
        ...     [[subsection]]
        ...     more_options = False
        ...     # end of file'''.splitlines()
        >>> b = '''# File is user.ini
        ...     [section1]
        ...     option1 = False
        ...     # end of file'''.splitlines()
        >>> c1 = ConfigObj(b)
        >>> c2 = ConfigObj(a)
        >>> c2.merge(c1)
        >>> c2
        ConfigObj({'section1': {'option1': 'False', 'subsection': {'more_options': 'False'}}})
        """
        ...
    
    def rename(self, oldkey, newkey): # -> None:
        """
        Change a keyname to another, without changing position in sequence.

        Implemented so that transformations can be made on keys,
        as well as on values. (used by encode and decode)

        Also renames comments.
        """
        ...
    
    def walk(self, function, raise_errors=..., call_on_sections=..., **keywargs): # -> dict[Unknown, Unknown]:
        """
        Walk every member and call a function on the keyword and value.

        Return a dictionary of the return values

        If the function raises an exception, raise the errror
        unless ``raise_errors=False``, in which case set the return value to
        ``False``.

        Any unrecognized keyword arguments you pass to walk, will be pased on
        to the function you pass in.

        Note: if ``call_on_sections`` is ``True`` then - on encountering a
        subsection, *first* the function is called for the *whole* subsection,
        and then recurses into it's members. This means your function must be
        able to handle strings, dictionaries and lists. This allows you
        to change the key of subsections as well as for ordinary members. The
        return value when called on the whole subsection has to be discarded.

        See  the encode and decode methods for examples, including functions.

        .. admonition:: caution

            You can use ``walk`` to transform the names of members of a section
            but you mustn't add or delete members.

        >>> config = '''[XXXXsection]
        ... XXXXkey = XXXXvalue'''.splitlines()
        >>> cfg = ConfigObj(config)
        >>> cfg
        ConfigObj({'XXXXsection': {'XXXXkey': 'XXXXvalue'}})
        >>> def transform(section, key):
        ...     val = section[key]
        ...     newkey = key.replace('XXXX', 'CLIENT1')
        ...     section.rename(key, newkey)
        ...     if isinstance(val, (tuple, list, dict)):
        ...         pass
        ...     else:
        ...         val = val.replace('XXXX', 'CLIENT1')
        ...         section[newkey] = val
        >>> cfg.walk(transform, call_on_sections=True)
        {'CLIENT1section': {'CLIENT1key': None}}
        >>> cfg
        ConfigObj({'CLIENT1section': {'CLIENT1key': 'CLIENT1value'}})
        """
        ...
    
    def as_bool(self, key): # -> bool:
        """
        Accepts a key as input. The corresponding value must be a string or
        the objects (``True`` or 1) or (``False`` or 0). We allow 0 and 1 to
        retain compatibility with Python 2.2.

        If the string is one of  ``True``, ``On``, ``Yes``, or ``1`` it returns
        ``True``.

        If the string is one of  ``False``, ``Off``, ``No``, or ``0`` it returns
        ``False``.

        ``as_bool`` is not case sensitive.

        Any other input will raise a ``ValueError``.

        >>> a = ConfigObj()
        >>> a['a'] = 'fish'
        >>> a.as_bool('a')
        Traceback (most recent call last):
        ValueError: Value "fish" is neither True nor False
        >>> a['b'] = 'True'
        >>> a.as_bool('b')
        1
        >>> a['b'] = 'off'
        >>> a.as_bool('b')
        0
        """
        ...
    
    def as_int(self, key): # -> int:
        """
        A convenience method which coerces the specified value to an integer.

        If the value is an invalid literal for ``int``, a ``ValueError`` will
        be raised.

        >>> a = ConfigObj()
        >>> a['a'] = 'fish'
        >>> a.as_int('a')
        Traceback (most recent call last):
        ValueError: invalid literal for int() with base 10: 'fish'
        >>> a['b'] = '1'
        >>> a.as_int('b')
        1
        >>> a['b'] = '3.2'
        >>> a.as_int('b')
        Traceback (most recent call last):
        ValueError: invalid literal for int() with base 10: '3.2'
        """
        ...
    
    def as_float(self, key): # -> float:
        """
        A convenience method which coerces the specified value to a float.

        If the value is an invalid literal for ``float``, a ``ValueError`` will
        be raised.

        >>> a = ConfigObj()
        >>> a['a'] = 'fish'
        >>> a.as_float('a')  #doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
        ValueError: invalid literal for float(): fish
        >>> a['b'] = '1'
        >>> a.as_float('b')
        1.0
        >>> a['b'] = '3.2'
        >>> a.as_float('b')  #doctest: +ELLIPSIS
        3.2...
        """
        ...
    
    def as_list(self, key): # -> list[str | Unknown]:
        """
        A convenience method which fetches the specified value, guaranteeing
        that it is a list.

        >>> a = ConfigObj()
        >>> a['a'] = 1
        >>> a.as_list('a')
        [1]
        >>> a['a'] = (1,)
        >>> a.as_list('a')
        [1]
        >>> a['a'] = [1]
        >>> a.as_list('a')
        [1]
        """
        ...
    
    def restore_default(self, key):
        """
        Restore (and return) default value for the specified key.

        This method will only work for a ConfigObj that was created
        with a configspec and has been validated.

        If there is no default value for this key, ``KeyError`` is raised.
        """
        ...
    
    def restore_defaults(self): # -> None:
        """
        Recursively restore default values to all members
        that have them.

        This method will only work for a ConfigObj that was created
        with a configspec and has been validated.

        It doesn't delete or modify entries without default values.
        """
        ...
    


class ConfigObj(Section):
    """An object to read, create, and write config files."""
    _keyword = ...
    _sectionmarker = ...
    _valueexp = ...
    _listvalueexp = ...
    _nolistvalue = ...
    _single_line_single = ...
    _single_line_double = ...
    _multi_line_single = ...
    _multi_line_double = ...
    _triple_quote = ...
    _bools = ...
    def __init__(self, infile=..., options=..., configspec=..., encoding=..., interpolation=..., raise_errors=..., list_values=..., create_empty=..., file_error=..., stringify=..., indent_type=..., default_encoding=..., unrepr=..., write_empty_values=..., _inspec=...) -> None:
        """
        Parse a config file or create a config file object.

        ``ConfigObj(infile=None, configspec=None, encoding=None,
                    interpolation=True, raise_errors=False, list_values=True,
                    create_empty=False, file_error=False, stringify=True,
                    indent_type=None, default_encoding=None, unrepr=False,
                    write_empty_values=False, _inspec=False)``
        """
        ...
    
    def __repr__(self): # -> str:
        ...
    
    def write(self, outfile=..., section=...):
        """
        Write the current ConfigObj as a file

        tekNico: FIXME: use StringIO instead of real files

        >>> filename = a.filename
        >>> a.filename = 'test.ini'
        >>> a.write()
        >>> a.filename = filename
        >>> a == ConfigObj('test.ini', raise_errors=True)
        1
        >>> import os
        >>> os.remove('test.ini')
        """
        ...
    
    def validate(self, validator, preserve_errors=..., copy=..., section=...):
        """
        Test the ConfigObj against a configspec.

        It uses the ``validator`` object from *validate.py*.

        To run ``validate`` on the current ConfigObj, call: ::

            test = config.validate(validator)

        (Normally having previously passed in the configspec when the ConfigObj
        was created - you can dynamically assign a dictionary of checks to the
        ``configspec`` attribute of a section though).

        It returns ``True`` if everything passes, or a dictionary of
        pass/fails (True/False). If every member of a subsection passes, it
        will just have the value ``True``. (It also returns ``False`` if all
        members fail).

        In addition, it converts the values from strings to their native
        types if their checks pass (and ``stringify`` is set).

        If ``preserve_errors`` is ``True`` (``False`` is default) then instead
        of a marking a fail with a ``False``, it will preserve the actual
        exception object. This can contain info about the reason for failure.
        For example the ``VdtValueTooSmallError`` indicates that the value
        supplied was too small. If a value (or section) is missing it will
        still be marked as ``False``.

        You must have the validate module to use ``preserve_errors=True``.

        You can then use the ``flatten_errors`` function to turn your nested
        results dictionary into a flattened list of failures - useful for
        displaying meaningful error messages.
        """
        ...
    
    def reset(self): # -> None:
        """Clear ConfigObj instance and restore to 'freshly created' state."""
        ...
    
    def reload(self): # -> None:
        """
        Reload a ConfigObj from file.

        This method raises a ``ReloadError`` if the ConfigObj doesn't have
        a filename attribute pointing to a file.
        """
        ...
    


class SimpleVal:
    """
    A simple validator.
    Can be used to check that all members expected are present.

    To use it, provide a configspec with all your members in (the value given
    will be ignored). Pass an instance of ``SimpleVal`` to the ``validate``
    method of your ``ConfigObj``. ``validate`` will return ``True`` if all
    members are present, or a dictionary with True/False meaning
    present/missing. (Whole missing sections will be replaced with ``False``)
    """
    def __init__(self) -> None:
        ...
    
    def check(self, check, member, missing=...):
        """A dummy check method, always returns the value unchanged."""
        ...
    


def flatten_errors(cfg, res, levels=..., results=...): # -> list[Unknown]:
    """
    An example function that will turn a nested dictionary of results
    (as returned by ``ConfigObj.validate``) into a flat list.

    ``cfg`` is the ConfigObj instance being checked, ``res`` is the results
    dictionary returned by ``validate``.

    (This is a recursive function, so you shouldn't use the ``levels`` or
    ``results`` arguments - they are used by the function.)

    Returns a list of keys that failed. Each member of the list is a tuple::

        ([list of sections...], key, result)

    If ``validate`` was called with ``preserve_errors=False`` (the default)
    then ``result`` will always be ``False``.

    *list of sections* is a flattened list of sections that the key was found
    in.

    If the section was missing (or a section was expected and a scalar provided
    - or vice-versa) then key will be ``None``.

    If the value (or section) was missing then ``result`` will be ``False``.

    If ``validate`` was called with ``preserve_errors=True`` and a value
    was present, but failed the check, then ``result`` will be the exception
    object returned. You can use this as a string that describes the failure.

    For example *The value "3" is of the wrong type*.
    """
    ...

def get_extra_values(conf, _prepend=...): # -> list[Unknown]:
    """
    Find all the values and sections not in the configspec from a validated
    ConfigObj.

    ``get_extra_values`` returns a list of tuples where each tuple represents
    either an extra section, or an extra value.

    The tuples contain two values, a tuple representing the section the value
    is in and the name of the extra values. For extra values in the top level
    section the first member will be an empty tuple. For values in the 'foo'
    section the first member will be ``('foo',)``. For members in the 'bar'
    subsection of the 'foo' section the first member will be ``('foo', 'bar')``.

    NOTE: If you call ``get_extra_values`` on a ConfigObj instance that hasn't
    been validated it will return an empty list.
    """
    ...

