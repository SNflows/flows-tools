"""
This type stub file was generated by pyright.
"""

import sys
from astropy.utils.decorators import deprecated_renamed_argument

"""Functions related to Python runtime introspection."""
__all__ = ['resolve_name', 'minversion', 'find_current_module', 'isinstancemethod']
__doctest_skip__ = ...
if sys.version_info[: 2] >= (3, 10):
    ...
else:
    def packages_distributions(): # -> dict[Unknown, list[Unknown]]:
        """
        Return a mapping of top-level packages to their distributions.
        Note: copied from https://github.com/python/importlib_metadata/pull/287
        """
        ...
    
def resolve_name(name, *additional_parts): # -> Any | Literal['']:
    """Resolve a name like ``module.object`` to an object and return it.

    This ends up working like ``from module import object`` but is easier
    to deal with than the `__import__` builtin and supports digging into
    submodules.

    Parameters
    ----------

    name : `str`
        A dotted path to a Python object--that is, the name of a function,
        class, or other object in a module with the full path to that module,
        including parent modules, separated by dots.  Also known as the fully
        qualified name of the object.

    additional_parts : iterable, optional
        If more than one positional arguments are given, those arguments are
        automatically dotted together with ``name``.

    Examples
    --------

    >>> resolve_name('astropy.utils.introspection.resolve_name')
    <function resolve_name at 0x...>
    >>> resolve_name('astropy', 'utils', 'introspection', 'resolve_name')
    <function resolve_name at 0x...>

    Raises
    ------
    `ImportError`
        If the module or named object is not found.
    """
    ...

@deprecated_renamed_argument('version_path', None, '5.0')
def minversion(module, version, inclusive=..., version_path=...): # -> bool:
    """
    Returns `True` if the specified Python module satisfies a minimum version
    requirement, and `False` if not.

    .. deprecated::
        ``version_path`` is not used anymore and is deprecated in
        ``astropy`` 5.0.

    Parameters
    ----------
    module : module or `str`
        An imported module of which to check the version, or the name of
        that module (in which case an import of that module is attempted--
        if this fails `False` is returned).

    version : `str`
        The version as a string that this module must have at a minimum (e.g.
        ``'0.12'``).

    inclusive : `bool`
        The specified version meets the requirement inclusively (i.e. ``>=``)
        as opposed to strictly greater than (default: `True`).

    Examples
    --------

    >>> import astropy
    >>> minversion(astropy, '0.4.4')
    True
    """
    ...

def find_current_module(depth=..., finddiff=...): # -> ModuleType | None:
    """
    Determines the module/package from which this function is called.

    This function has two modes, determined by the ``finddiff`` option. it
    will either simply go the requested number of frames up the call
    stack (if ``finddiff`` is False), or it will go up the call stack until
    it reaches a module that is *not* in a specified set.

    Parameters
    ----------
    depth : int
        Specifies how far back to go in the call stack (0-indexed, so that
        passing in 0 gives back `astropy.utils.misc`).
    finddiff : bool or list
        If False, the returned ``mod`` will just be ``depth`` frames up from
        the current frame. Otherwise, the function will start at a frame
        ``depth`` up from current, and continue up the call stack to the
        first module that is *different* from those in the provided list.
        In this case, ``finddiff`` can be a list of modules or modules
        names. Alternatively, it can be True, which will use the module
        ``depth`` call stack frames up as the module the returned module
        most be different from.

    Returns
    -------
    mod : module or None
        The module object or None if the package cannot be found. The name of
        the module is available as the ``__name__`` attribute of the returned
        object (if it isn't None).

    Raises
    ------
    ValueError
        If ``finddiff`` is a list with an invalid entry.

    Examples
    --------
    The examples below assume that there are two modules in a package named
    ``pkg``. ``mod1.py``::

        def find1():
            from astropy.utils import find_current_module
            print find_current_module(1).__name__
        def find2():
            from astropy.utils import find_current_module
            cmod = find_current_module(2)
            if cmod is None:
                print 'None'
            else:
                print cmod.__name__
        def find_diff():
            from astropy.utils import find_current_module
            print find_current_module(0,True).__name__

    ``mod2.py``::

        def find():
            from .mod1 import find2
            find2()

    With these modules in place, the following occurs::

        >>> from pkg import mod1, mod2
        >>> from astropy.utils import find_current_module
        >>> mod1.find1()
        pkg.mod1
        >>> mod1.find2()
        None
        >>> mod2.find()
        pkg.mod2
        >>> find_current_module(0)
        <module 'astropy.utils.misc' from 'astropy/utils/misc.py'>
        >>> mod1.find_diff()
        pkg.mod1

    """
    ...

def find_mod_objs(modname, onlylocals=...): # -> tuple[list[Any | Unknown | str], list[Unknown], list[Any]]:
    """ Returns all the public attributes of a module referenced by name.

    .. note::
        The returned list *not* include subpackages or modules of
        ``modname``, nor does it include private attributes (those that
        begin with '_' or are not in `__all__`).

    Parameters
    ----------
    modname : str
        The name of the module to search.
    onlylocals : bool or list of str
        If `True`, only attributes that are either members of ``modname`` OR
        one of its modules or subpackages will be included. If it is a list
        of strings, those specify the possible packages that will be
        considered "local".

    Returns
    -------
    localnames : list of str
        A list of the names of the attributes as they are named in the
        module ``modname`` .
    fqnames : list of str
        A list of the full qualified names of the attributes (e.g.,
        ``astropy.utils.introspection.find_mod_objs``). For attributes that are
        simple variables, this is based on the local name, but for functions or
        classes it can be different if they are actually defined elsewhere and
        just referenced in ``modname``.
    objs : list of objects
        A list of the actual attributes themselves (in the same order as
        the other arguments)

    """
    ...

def isinstancemethod(cls, obj): # -> bool:
    """
    Returns `True` if the given object is an instance method of the class
    it is defined on (as opposed to a `staticmethod` or a `classmethod`).

    This requires both the class the object is a member of as well as the
    object itself in order to make this determination.

    Parameters
    ----------
    cls : `type`
        The class on which this method was defined.
    obj : `object`
        A member of the provided class (the membership is not checked directly,
        but this function will always return `False` if the given object is not
        a member of the given class).

    Examples
    --------
    >>> class MetaClass(type):
    ...     def a_classmethod(cls): pass
    ...
    >>> class MyClass(metaclass=MetaClass):
    ...     def an_instancemethod(self): pass
    ...
    ...     @classmethod
    ...     def another_classmethod(cls): pass
    ...
    ...     @staticmethod
    ...     def a_staticmethod(): pass
    ...
    >>> isinstancemethod(MyClass, MyClass.a_classmethod)
    False
    >>> isinstancemethod(MyClass, MyClass.another_classmethod)
    False
    >>> isinstancemethod(MyClass, MyClass.a_staticmethod)
    False
    >>> isinstancemethod(MyClass, MyClass.an_instancemethod)
    True
    """
    ...

