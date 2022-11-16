"""
This type stub file was generated by pyright.
"""

import abc
import contextlib

__all__ = ['IORegistryError']
class IORegistryError(Exception):
    """Custom error for registry clashes.
    """
    ...


class _UnifiedIORegistryBase(metaclass=abc.ABCMeta):
    """Base class for registries in Astropy's Unified IO.

    This base class provides identification functions and miscellaneous
    utilities. For an example how to build a registry subclass we suggest
    :class:`~astropy.io.registry.UnifiedInputRegistry`, which enables
    read-only registries. These higher-level subclasses will probably serve
    better as a baseclass, for instance
    :class:`~astropy.io.registry.UnifiedIORegistry` subclasses both
    :class:`~astropy.io.registry.UnifiedInputRegistry` and
    :class:`~astropy.io.registry.UnifiedOutputRegistry` to enable both
    reading from and writing to files.

    .. versionadded:: 5.0

    """
    def __init__(self) -> None:
        ...
    
    @property
    def available_registries(self): # -> dict_keys[Unknown, Unknown]:
        """Available registries.

        Returns
        -------
        ``dict_keys``
        """
        ...
    
    def get_formats(self, data_class=..., filter_on=...): # -> Table:
        """
        Get the list of registered formats as a `~astropy.table.Table`.

        Parameters
        ----------
        data_class : class or None, optional
            Filter readers/writer to match data class (default = all classes).
        filter_on : str or None, optional
            Which registry to show. E.g. "identify"
            If None search for both.  Default is None.

        Returns
        -------
        format_table : :class:`~astropy.table.Table`
            Table of available I/O formats.

        Raises
        ------
        ValueError
            If ``filter_on`` is not None nor a registry name.
        """
        ...
    
    @contextlib.contextmanager
    def delay_doc_updates(self, cls): # -> Generator[None, None, None]:
        """Contextmanager to disable documentation updates when registering
        reader and writer. The documentation is only built once when the
        contextmanager exits.

        .. versionadded:: 1.3

        Parameters
        ----------
        cls : class
            Class for which the documentation updates should be delayed.

        Notes
        -----
        Registering multiple readers and writers can cause significant overhead
        because the documentation of the corresponding ``read`` and ``write``
        methods are build every time.

        Examples
        --------
        see for example the source code of ``astropy.table.__init__``.
        """
        ...
    
    def register_identifier(self, data_format, data_class, identifier, force=...): # -> None:
        """
        Associate an identifier function with a specific data type.

        Parameters
        ----------
        data_format : str
            The data format identifier. This is the string that is used to
            specify the data type when reading/writing.
        data_class : class
            The class of the object that can be written.
        identifier : function
            A function that checks the argument specified to `read` or `write` to
            determine whether the input can be interpreted as a table of type
            ``data_format``. This function should take the following arguments:

               - ``origin``: A string ``"read"`` or ``"write"`` identifying whether
                 the file is to be opened for reading or writing.
               - ``path``: The path to the file.
               - ``fileobj``: An open file object to read the file's contents, or
                 `None` if the file could not be opened.
               - ``*args``: Positional arguments for the `read` or `write`
                 function.
               - ``**kwargs``: Keyword arguments for the `read` or `write`
                 function.

            One or both of ``path`` or ``fileobj`` may be `None`.  If they are
            both `None`, the identifier will need to work from ``args[0]``.

            The function should return True if the input can be identified
            as being of format ``data_format``, and False otherwise.
        force : bool, optional
            Whether to override any existing function if already present.
            Default is ``False``.

        Examples
        --------
        To set the identifier based on extensions, for formats that take a
        filename as a first argument, you can do for example

        .. code-block:: python

            from astropy.io.registry import register_identifier
            from astropy.table import Table
            def my_identifier(*args, **kwargs):
                return isinstance(args[0], str) and args[0].endswith('.tbl')
            register_identifier('ipac', Table, my_identifier)
            unregister_identifier('ipac', Table)
        """
        ...
    
    def unregister_identifier(self, data_format, data_class): # -> None:
        """
        Unregister an identifier function

        Parameters
        ----------
        data_format : str
            The data format identifier.
        data_class : class
            The class of the object that can be read/written.
        """
        ...
    
    def identify_format(self, origin, data_class_required, path, fileobj, args, kwargs): # -> list[Unknown]:
        """Loop through identifiers to see which formats match.

        Parameters
        ----------
        origin : str
            A string ``"read`` or ``"write"`` identifying whether the file is to be
            opened for reading or writing.
        data_class_required : object
            The specified class for the result of `read` or the class that is to be
            written.
        path : str or path-like or None
            The path to the file or None.
        fileobj : file-like or None.
            An open file object to read the file's contents, or ``None`` if the
            file could not be opened.
        args : sequence
            Positional arguments for the `read` or `write` function. Note that
            these must be provided as sequence.
        kwargs : dict-like
            Keyword arguments for the `read` or `write` function. Note that this
            parameter must be `dict`-like.

        Returns
        -------
        valid_formats : list
            List of matching formats.
        """
        ...
    

