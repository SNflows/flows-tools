"""
This type stub file was generated by pyright.
"""

from astropy.utils.compat.optional_deps import HAS_BZ2

if HAS_BZ2:
    ...
IO_FITS_MODES = ...
FILE_MODES = ...
TEXT_RE = ...
MEMMAP_MODES = ...
GZIP_MAGIC = ...
PKZIP_MAGIC = ...
BZIP2_MAGIC = ...
class _File:
    """
    Represents a FITS file on disk (or in some other file-like object).
    """
    def __init__(self, fileobj=..., mode=..., memmap=..., overwrite=..., cache=...) -> None:
        ...
    
    def __repr__(self): # -> str:
        ...
    
    def __enter__(self): # -> Self@_File:
        ...
    
    def __exit__(self, type, value, traceback): # -> None:
        ...
    
    def readable(self): # -> bool:
        ...
    
    def read(self, size=...): # -> bytes | Any | Literal['']:
        ...
    
    def readarray(self, size=..., offset=..., dtype=..., shape=...): # -> ndarray[Unknown, Unknown] | NDArray[float64]:
        """
        Similar to file.read(), but returns the contents of the underlying
        file as a numpy array (or mmap'd array if memmap=True) rather than a
        string.

        Usually it's best not to use the `size` argument with this method, but
        it's provided for compatibility.
        """
        ...
    
    def writable(self): # -> bool:
        ...
    
    def write(self, string): # -> None:
        ...
    
    def writearray(self, array): # -> None:
        """
        Similar to file.write(), but writes a numpy array instead of a string.

        Also like file.write(), a flush() or close() may be needed before
        the file on disk reflects the data written.
        """
        ...
    
    def flush(self): # -> None:
        ...
    
    def seek(self, offset, whence=...): # -> None:
        ...
    
    def tell(self): # -> int:
        ...
    
    def truncate(self, size=...): # -> None:
        ...
    
    def close(self): # -> None:
        """
        Close the 'physical' FITS file.
        """
        ...
    


