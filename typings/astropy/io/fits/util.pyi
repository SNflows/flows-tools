"""
This type stub file was generated by pyright.
"""

path_like = ...
cmp = ...
all_integer_types = ...
class NotifierMixin:
    """
    Mixin class that provides services by which objects can register
    listeners to changes on that object.

    All methods provided by this class are underscored, since this is intended
    for internal use to communicate between classes in a generic way, and is
    not machinery that should be exposed to users of the classes involved.

    Use the ``_add_listener`` method to register a listener on an instance of
    the notifier.  This registers the listener with a weak reference, so if
    no other references to the listener exist it is automatically dropped from
    the list and does not need to be manually removed.

    Call the ``_notify`` method on the notifier to update all listeners
    upon changes.  ``_notify('change_type', *args, **kwargs)`` results
    in calling ``listener._update_change_type(*args, **kwargs)`` on all
    listeners subscribed to that notifier.

    If a particular listener does not have the appropriate update method
    it is ignored.

    Examples
    --------

    >>> class Widget(NotifierMixin):
    ...     state = 1
    ...     def __init__(self, name):
    ...         self.name = name
    ...     def update_state(self):
    ...         self.state += 1
    ...         self._notify('widget_state_changed', self)
    ...
    >>> class WidgetListener:
    ...     def _update_widget_state_changed(self, widget):
    ...         print('Widget {0} changed state to {1}'.format(
    ...             widget.name, widget.state))
    ...
    >>> widget = Widget('fred')
    >>> listener = WidgetListener()
    >>> widget._add_listener(listener)
    >>> widget.update_state()
    Widget fred changed state to 2
    """
    _listeners = ...
    def __getstate__(self): # -> dict[str, Any]:
        """
        Exclude listeners when saving the listener's state, since they may be
        ephemeral.
        """
        ...
    


def first(iterable):
    """
    Returns the first item returned by iterating over an iterable object.

    Example:

    >>> a = [1, 2, 3]
    >>> first(a)
    1
    """
    ...

def itersubclasses(cls, _seen=...): # -> Generator[Unknown, None, None]:
    """
    Generator over all subclasses of a given class, in depth first order.

    >>> class A: pass
    >>> class B(A): pass
    >>> class C(A): pass
    >>> class D(B,C): pass
    >>> class E(D): pass
    >>>
    >>> for cls in itersubclasses(A):
    ...     print(cls.__name__)
    B
    D
    E
    C
    >>> # get ALL classes currently defined
    >>> [cls.__name__ for cls in itersubclasses(object)]
    [...'tuple', ...'type', ...]

    From http://code.activestate.com/recipes/576949/
    """
    ...

def ignore_sigint(func): # -> (*args: Unknown, **kwargs: Unknown) -> None:
    """
    This decorator registers a custom SIGINT handler to catch and ignore SIGINT
    until the wrapped function is completed.
    """
    ...

def pairwise(iterable): # -> zip[tuple[Unknown, Unknown]]:
    """Return the items of an iterable paired with its next item.

    Ex: s -> (s0,s1), (s1,s2), (s2,s3), ....
    """
    ...

def encode_ascii(s): # -> bytes | NDArray[Any] | ndarray[Unknown, Unknown]:
    ...

def decode_ascii(s): # -> str | NDArray[Any] | ndarray[Unknown, Unknown]:
    ...

def isreadable(f): # -> bool:
    """
    Returns True if the file-like object can be read from.  This is a common-
    sense approximation of io.IOBase.readable.
    """
    ...

def iswritable(f): # -> bool:
    """
    Returns True if the file-like object can be written to.  This is a common-
    sense approximation of io.IOBase.writable.
    """
    ...

def isfile(f): # -> bool:
    """
    Returns True if the given object represents an OS-level file (that is,
    ``isinstance(f, file)``).

    On Python 3 this also returns True if the given object is higher level
    wrapper on top of a FileIO object, such as a TextIOWrapper.
    """
    ...

def fileobj_name(f): # -> str | bytes:
    """
    Returns the 'name' of file-like object *f*, if it has anything that could be
    called its name.  Otherwise f's class or type is returned.  If f is a
    string f itself is returned.
    """
    ...

def fileobj_closed(f): # -> bool:
    """
    Returns True if the given file-like object is closed or if *f* is a string
    (and assumed to be a pathname).

    Returns False for all other types of objects, under the assumption that
    they are file-like objects with no sense of a 'closed' state.
    """
    ...

def fileobj_mode(f): # -> Literal['rb', 'wb'] | None:
    """
    Returns the 'mode' string of a file-like object if such a thing exists.
    Otherwise returns None.
    """
    ...

def fileobj_is_binary(f): # -> bool:
    """
    Returns True if the give file or file-like object has a file open in binary
    mode.  When in doubt, returns True by default.
    """
    ...

def translate(s, table, deletechars):
    ...

def fill(text, width, **kwargs): # -> str:
    """
    Like :func:`textwrap.wrap` but preserves existing paragraphs which
    :func:`textwrap.wrap` does not otherwise handle well.  Also handles section
    headers.
    """
    ...

CHUNKED_FROMFILE = ...
_OSX_WRITE_LIMIT = (2 ** 32) - 1
_WIN_WRITE_LIMIT = (2 ** 31) - 1
def get_testdata_filepath(filename): # -> str | LiteralString:
    """
    Return a string representing the path to the file requested from the
    io.fits test data set.

    .. versionadded:: 2.0.3

    Parameters
    ----------
    filename : str
        The filename of the test data file.

    Returns
    -------
    filepath : str
        The path to the requested file.
    """
    ...

