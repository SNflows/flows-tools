"""
This type stub file was generated by pyright.
"""

from astropy import coordinates as coord
from astropy.utils import deprecated

"""
Common functions and classes that are required by all query classes.
"""
def ICRSCoordGenerator(*args, **kwargs): # -> SkyCoord:
    ...

def GalacticCoordGenerator(*args, **kwargs): # -> SkyCoord:
    ...

def FK5CoordGenerator(*args, **kwargs): # -> SkyCoord:
    ...

def FK4CoordGenerator(*args, **kwargs): # -> SkyCoord:
    ...

ICRSCoord = coord.SkyCoord
CoordClasses = ...
__all__ = ['send_request', 'parse_coordinates', 'TableList', 'suppress_vo_warnings', 'validate_email', 'ASTROPY_LT_4_1', 'ASTROPY_LT_4_3', 'ASTROPY_LT_5_0']
ASTROPY_LT_4_1 = ...
ASTROPY_LT_4_3 = ...
ASTROPY_LT_5_0 = ...
@deprecated('0.4.4', alternative='astroquery.query.BaseQuery._request')
def send_request(url, data, timeout, request_type=..., headers=..., **kwargs): # -> Response:
    """
    A utility function that post HTTP requests to remote server
    and returns the HTTP response.

    Parameters
    ----------
    url : str
        The URL of the remote server
    data : dict
        A dictionary representing the payload to be posted via the HTTP request
    timeout : int, quantity_like
        Time limit for establishing successful connection with remote server
    request_type : str
        options are 'POST' (default) and 'GET'. Determines whether to perform
        an HTTP POST or an HTTP GET request
    headers : dict
        POST or GET headers.  user-agent will be set to
        astropy:astroquery.version

    Returns
    -------
    response : `requests.Response`
        Response object returned by the remote server
    """
    ...

def radius_to_unit(radius, unit=...): # -> Any:
    """
    Helper function: Parse a radius, then return its value in degrees

    Parameters
    ----------
    radius : str or `~astropy.units.Quantity`
        The radius of a region

    Returns
    -------
    Floating point scalar value of radius in degrees
    """
    ...

def parse_coordinates(coordinates): # -> SkyCoord | BaseCoordinateFrame:
    """
    Takes a string or astropy.coordinates object. Checks if the
    string is parsable as an `astropy.coordinates`
    object or is a name that is resolvable. Otherwise asserts
    that the argument is an astropy.coordinates object.

    Parameters
    ----------
    coordinates : str or `astropy.coordinates` object
        Astronomical coordinate

    Returns
    -------
    coordinates : a subclass of `astropy.coordinates.BaseCoordinateFrame`


    Raises
    ------
    astropy.units.UnitsError
    TypeError
    """
    ...

def coord_to_radec(coordinate): # -> tuple[Unknown, Unknown]:
    """
    Wrapper to turn any astropy coordinate into FK5 RA in Hours and FK5 Dec in
    degrees

    This is a hack / temporary wrapper to deal with the unstable astropy API
    (it may be wise to remove this hack since it's not clear that the old
    coordinate API can even do transforms)
    """
    ...

class TableList(list):
    """
    A class that inherits from `list` but included some pretty printing methods
    for an OrderedDict of `astropy.table.Table` objects.

    HINT: To access the tables by # instead of by table ID:
    >>> t = TableList([('a',1),('b',2)])
    >>> t[1]
    2
    >>> t['b']
    2
    """
    def __init__(self, inp) -> None:
        ...
    
    def __getitem__(self, key):
        ...
    
    def __setitem__(self, value):
        ...
    
    def __getslice__(self, slice):
        ...
    
    def keys(self): # -> list[Unknown]:
        ...
    
    def values(self): # -> list[Unknown]:
        ...
    
    def __repr__(self): # -> str:
        """
        Overrides the `OrderedDict.__repr__` method to return a simple summary
        of the `TableList` object.
        """
        ...
    
    def format_table_list(self): # -> str:
        """
        Prints the names of all `astropy.table.Table` objects, with their
        respective number of row and columns, contained in the
        `TableList` instance.
        """
        ...
    
    def print_table_list(self): # -> None:
        ...
    
    def pprint(self, **kwargs): # -> None:
        """ Helper function to make API more similar to astropy.Tables """
        ...
    


def suppress_vo_warnings(): # -> None:
    """
    Suppresses all warnings of the class
    `astropy.io.votable.exceptions.VOWarning`.
    """
    ...

def validate_email(email): # -> bool:
    """
    E-mail address validation.  Uses validate_email if available, else a simple
    regex that will let through some invalid e-mails but will catch the most
    common violators.
    """
    ...

class FileContainer:
    """
    A File Object container, meant to offer lazy access to downloaded FITS
    files.
    """
    def __init__(self, target, **kwargs) -> None:
        ...
    
    def get_fits(self):
        """
        Assuming the contained file is a FITS file, read it
        and return the file parsed as FITS HDUList
        """
        ...
    
    def save_fits(self, savepath, link_cache=...): # -> None:
        """
        Save a FITS file to savepath

        Parameters
        ----------
        savepath : str
            The full path to a FITS filename, e.g. "file.fits", or
            "/path/to/file.fits".
        link_cache : 'hard', 'sym', or False
            Try to create a hard or symbolic link to the astropy cached file?
            If the system is unable to create a hardlink, the file will be
            copied to the target location.
        """
        ...
    
    def get_string(self):
        """
        Download the file as a string
        """
        ...
    
    def get_stringio(self): # -> BytesIO | StringIO:
        """
        Return the file as an io.StringIO object
        """
        ...
    
    def __repr__(self): # -> str:
        ...
    


def get_readable_fileobj(*args, **kwargs): # -> _GeneratorContextManager[Unknown]:
    """
    Overload astropy's get_readable_fileobj so that we can safely monkeypatch
    it in astroquery without affecting astropy core functionality
    """
    ...

def parse_votable(content):
    """
    Parse a votable in string format
    """
    ...

