"""
This type stub file was generated by pyright.
"""

from .decorators import classproperty

"""
Utilities for console input and output.
"""
__all__ = ['isatty', 'color_print', 'human_time', 'human_file_size', 'ProgressBar', 'Spinner', 'print_code_line', 'ProgressBarOrSpinner', 'terminal_size']
_DEFAULT_ENCODING = ...
class _IPython:
    """Singleton class given access to IPython streams, etc."""
    @classproperty
    def get_ipython(cls):
        ...
    
    @classproperty
    def OutStream(cls): # -> None:
        ...
    
    @classproperty
    def ipyio(cls): # -> None:
        ...
    
    @classmethod
    def get_stream(cls, stream): # -> Any:
        ...
    


def isatty(file): # -> bool:
    """
    Returns `True` if ``file`` is a tty.

    Most built-in Python file-like objects have an `isatty` member,
    but some user-defined types may not, so this assumes those are not
    ttys.
    """
    ...

def terminal_size(file=...): # -> tuple[Any, Any] | tuple[int, int] | tuple[ConfigItem | Unknown | int | bool | list[int | float | str | bool] | string_type, ConfigItem | Unknown | int | bool | list[int | float | str | bool] | string_type]:
    """
    Returns a tuple (height, width) containing the height and width of
    the terminal.

    This function will look for the width in height in multiple areas
    before falling back on the width and height in astropy's
    configuration.
    """
    ...

def color_print(*args, end=..., **kwargs): # -> None:
    """
    Prints colors and styles to the terminal uses ANSI escape
    sequences.

    ::

       color_print('This is the color ', 'default', 'GREEN', 'green')

    Parameters
    ----------
    positional args : str
        The positional arguments come in pairs (*msg*, *color*), where
        *msg* is the string to display and *color* is the color to
        display it in.

        *color* is an ANSI terminal color name.  Must be one of:
        black, red, green, brown, blue, magenta, cyan, lightgrey,
        default, darkgrey, lightred, lightgreen, yellow, lightblue,
        lightmagenta, lightcyan, white, or '' (the empty string).

    file : writable file-like, optional
        Where to write to.  Defaults to `sys.stdout`.  If file is not
        a tty (as determined by calling its `isatty` member, if one
        exists), no coloring will be included.

    end : str, optional
        The ending of the message.  Defaults to ``\\n``.  The end will
        be printed after resetting any color or font state.
    """
    ...

def strip_ansi_codes(s): # -> str:
    """
    Remove ANSI color codes from the string.
    """
    ...

def human_time(seconds): # -> str:
    """
    Returns a human-friendly time string that is always exactly 6
    characters long.

    Depending on the number of seconds given, can be one of::

        1w 3d
        2d 4h
        1h 5m
        1m 4s
          15s

    Will be in color if console coloring is turned on.

    Parameters
    ----------
    seconds : int
        The number of seconds to represent

    Returns
    -------
    time : str
        A human-friendly representation of the given number of seconds
        that is always exactly 6 characters.
    """
    ...

def human_file_size(size): # -> str:
    """
    Returns a human-friendly string representing a file size
    that is 2-4 characters long.

    For example, depending on the number of bytes given, can be one
    of::

        256b
        64k
        1.1G

    Parameters
    ----------
    size : int
        The size of the file (in bytes)

    Returns
    -------
    size : str
        A human-friendly representation of the size of the file
    """
    ...

class _mapfunc:
    """
    A function wrapper to support ProgressBar.map().
    """
    def __init__(self, func) -> None:
        ...
    
    def __call__(self, i_arg): # -> tuple[Unknown, Unknown]:
        ...
    


class ProgressBar:
    """
    A class to display a progress bar in the terminal.

    It is designed to be used either with the ``with`` statement::

        with ProgressBar(len(items)) as bar:
            for item in enumerate(items):
                bar.update()

    or as a generator::

        for item in ProgressBar(items):
            item.process()
    """
    def __init__(self, total_or_items, ipython_widget=..., file=...) -> None:
        """
        Parameters
        ----------
        total_or_items : int or sequence
            If an int, the number of increments in the process being
            tracked.  If a sequence, the items to iterate over.

        ipython_widget : bool, optional
            If `True`, the progress bar will display as an IPython
            notebook widget.

        file : writable file-like, optional
            The file to write the progress bar to.  Defaults to
            `sys.stdout`.  If ``file`` is not a tty (as determined by
            calling its `isatty` member, if any, or special case hacks
            to detect the IPython console), the progress bar will be
            completely silent.
        """
        ...
    
    def __enter__(self): # -> Self@ProgressBar:
        ...
    
    def __exit__(self, exc_type, exc_value, traceback): # -> None:
        ...
    
    def __iter__(self): # -> Self@ProgressBar:
        ...
    
    def __next__(self): # -> int:
        ...
    
    def update(self, value=...): # -> None:
        """
        Update progress bar via the console or notebook accordingly.
        """
        ...
    
    @classmethod
    def map(cls, function, items, multiprocess=..., file=..., step=..., ipython_widget=..., multiprocessing_start_method=...): # -> list[Unknown]:
        """Map function over items while displaying a progress bar with percentage complete.

        The map operation may run in arbitrary order on the items, but the results are
        returned in sequential order.

        ::

            def work(i):
                print(i)

            ProgressBar.map(work, range(50))

        Parameters
        ----------
        function : function
            Function to call for each step

        items : sequence
            Sequence where each element is a tuple of arguments to pass to
            *function*.

        multiprocess : bool, int, optional
            If `True`, use the `multiprocessing` module to distribute each task
            to a different processor core. If a number greater than 1, then use
            that number of cores.

        ipython_widget : bool, optional
            If `True`, the progress bar will display as an IPython
            notebook widget.

        file : writable file-like, optional
            The file to write the progress bar to.  Defaults to
            `sys.stdout`.  If ``file`` is not a tty (as determined by
            calling its `isatty` member, if any), the scrollbar will
            be completely silent.

        step : int, optional
            Update the progress bar at least every *step* steps (default: 100).
            If ``multiprocess`` is `True`, this will affect the size
            of the chunks of ``items`` that are submitted as separate tasks
            to the process pool.  A large step size may make the job
            complete faster if ``items`` is very long.

        multiprocessing_start_method : str, optional
            Useful primarily for testing; if in doubt leave it as the default.
            When using multiprocessing, certain anomalies occur when starting
            processes with the "spawn" method (the only option on Windows);
            other anomalies occur with the "fork" method (the default on
            Linux).
        """
        ...
    
    @classmethod
    def map_unordered(cls, function, items, multiprocess=..., file=..., step=..., ipython_widget=..., multiprocessing_start_method=...): # -> list[Unknown]:
        """Map function over items, reporting the progress.

        Does a `map` operation while displaying a progress bar with
        percentage complete. The map operation may run on arbitrary order
        on the items, and the results may be returned in arbitrary order.

        ::

            def work(i):
                print(i)

            ProgressBar.map(work, range(50))

        Parameters
        ----------
        function : function
            Function to call for each step

        items : sequence
            Sequence where each element is a tuple of arguments to pass to
            *function*.

        multiprocess : bool, int, optional
            If `True`, use the `multiprocessing` module to distribute each task
            to a different processor core. If a number greater than 1, then use
            that number of cores.

        ipython_widget : bool, optional
            If `True`, the progress bar will display as an IPython
            notebook widget.

        file : writable file-like, optional
            The file to write the progress bar to.  Defaults to
            `sys.stdout`.  If ``file`` is not a tty (as determined by
            calling its `isatty` member, if any), the scrollbar will
            be completely silent.

        step : int, optional
            Update the progress bar at least every *step* steps (default: 100).
            If ``multiprocess`` is `True`, this will affect the size
            of the chunks of ``items`` that are submitted as separate tasks
            to the process pool.  A large step size may make the job
            complete faster if ``items`` is very long.

        multiprocessing_start_method : str, optional
            Useful primarily for testing; if in doubt leave it as the default.
            When using multiprocessing, certain anomalies occur when starting
            processes with the "spawn" method (the only option on Windows);
            other anomalies occur with the "fork" method (the default on
            Linux).
        """
        ...
    


class Spinner:
    """
    A class to display a spinner in the terminal.

    It is designed to be used with the ``with`` statement::

        with Spinner("Reticulating splines", "green") as s:
            for item in enumerate(items):
                s.update()
    """
    _default_unicode_chars = ...
    _default_ascii_chars = ...
    def __init__(self, msg, color=..., file=..., step=..., chars=...) -> None:
        """
        Parameters
        ----------
        msg : str
            The message to print

        color : str, optional
            An ANSI terminal color name.  Must be one of: black, red,
            green, brown, blue, magenta, cyan, lightgrey, default,
            darkgrey, lightred, lightgreen, yellow, lightblue,
            lightmagenta, lightcyan, white.

        file : writable file-like, optional
            The file to write the spinner to.  Defaults to
            `sys.stdout`.  If ``file`` is not a tty (as determined by
            calling its `isatty` member, if any, or special case hacks
            to detect the IPython console), the spinner will be
            completely silent.

        step : int, optional
            Only update the spinner every *step* steps

        chars : str, optional
            The character sequence to use for the spinner
        """
        ...
    
    def __enter__(self): # -> Self@Spinner:
        ...
    
    def __exit__(self, exc_type, exc_value, traceback): # -> None:
        ...
    
    def __iter__(self): # -> Self@Spinner:
        ...
    
    def __next__(self): # -> None:
        ...
    
    def update(self, value=...): # -> None:
        """Update the spin wheel in the terminal.

        Parameters
        ----------
        value : int, optional
            Ignored (present just for compatibility with `ProgressBar.update`).

        """
        ...
    


class ProgressBarOrSpinner:
    """
    A class that displays either a `ProgressBar` or `Spinner`
    depending on whether the total size of the operation is
    known or not.

    It is designed to be used with the ``with`` statement::

        if file.has_length():
            length = file.get_length()
        else:
            length = None
        bytes_read = 0
        with ProgressBarOrSpinner(length) as bar:
            while file.read(blocksize):
                bytes_read += blocksize
                bar.update(bytes_read)
    """
    def __init__(self, total, msg, color=..., file=...) -> None:
        """
        Parameters
        ----------
        total : int or None
            If an int, the number of increments in the process being
            tracked and a `ProgressBar` is displayed.  If `None`, a
            `Spinner` is displayed.

        msg : str
            The message to display above the `ProgressBar` or
            alongside the `Spinner`.

        color : str, optional
            The color of ``msg``, if any.  Must be an ANSI terminal
            color name.  Must be one of: black, red, green, brown,
            blue, magenta, cyan, lightgrey, default, darkgrey,
            lightred, lightgreen, yellow, lightblue, lightmagenta,
            lightcyan, white.

        file : writable file-like, optional
            The file to write the to.  Defaults to `sys.stdout`.  If
            ``file`` is not a tty (as determined by calling its `isatty`
            member, if any), only ``msg`` will be displayed: the
            `ProgressBar` or `Spinner` will be silent.
        """
        ...
    
    def __enter__(self): # -> Self@ProgressBarOrSpinner:
        ...
    
    def __exit__(self, exc_type, exc_value, traceback): # -> None:
        ...
    
    def update(self, value): # -> None:
        """
        Update the progress bar to the given value (out of the total
        given to the constructor.
        """
        ...
    


def print_code_line(line, col=..., file=..., tabwidth=..., width=...): # -> None:
    """
    Prints a line of source code, highlighting a particular character
    position in the line.  Useful for displaying the context of error
    messages.

    If the line is more than ``width`` characters, the line is truncated
    accordingly and '…' characters are inserted at the front and/or
    end.

    It looks like this::

        there_is_a_syntax_error_here :
                                     ^

    Parameters
    ----------
    line : unicode
        The line of code to display

    col : int, optional
        The character in the line to highlight.  ``col`` must be less
        than ``len(line)``.

    file : writable file-like, optional
        Where to write to.  Defaults to `sys.stdout`.

    tabwidth : int, optional
        The number of spaces per tab (``'\\t'``) character.  Default
        is 8.  All tabs will be converted to spaces to ensure that the
        caret lines up with the correct column.

    width : int, optional
        The width of the display, beyond which the line will be
        truncated.  Defaults to 70 (this matches the default in the
        standard library's `textwrap` module).
    """
    ...

class Getch:
    """Get a single character from standard input without screen echo.

    Returns
    -------
    char : str (one character)
    """
    def __init__(self) -> None:
        ...
    
    def __call__(self): # -> str:
        ...
    


class _GetchUnix:
    def __init__(self) -> None:
        ...
    
    def __call__(self): # -> str:
        ...
    


class _GetchWindows:
    def __init__(self) -> None:
        ...
    
    def __call__(self):
        ...
    


class _GetchMacCarbon:
    """
    A function which returns the current ASCII key that is down;
    if no ASCII key is down, the null string is returned.  The
    page http://www.mactech.com/macintosh-c/chap02-1.html was
    very helpful in figuring out how to do this.
    """
    def __init__(self) -> None:
        ...
    
    def __call__(self): # -> str:
        ...
    


