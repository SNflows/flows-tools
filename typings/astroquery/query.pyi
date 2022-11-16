"""
This type stub file was generated by pyright.
"""

import abc

__all__ = ['BaseQuery', 'QueryWithLogin']
def to_cache(response, cache_file): # -> None:
    ...

class AstroQuery:
    def __init__(self, method, url, params=..., data=..., headers=..., files=..., timeout=..., json=...) -> None:
        ...
    
    @property
    def timeout(self):
        ...
    
    @timeout.setter
    def timeout(self, value): # -> None:
        ...
    
    def request(self, session, cache_location=..., stream=..., auth=..., verify=..., allow_redirects=..., json=...):
        ...
    
    def hash(self): # -> str:
        ...
    
    def request_file(self, cache_location): # -> str:
        ...
    
    def from_cache(self, cache_location): # -> Response | None:
        ...
    
    def remove_cache_file(self, cache_location): # -> None:
        """
        Remove the cache file - may be needed if a query fails during parsing
        (successful request, but failed return)
        """
        ...
    


class LoginABCMeta(abc.ABCMeta):
    """
    The goal of this metaclass is to copy the docstring and signature from
    ._login methods, implemented in subclasses, to a .login method that is
    visible by the users.

    It also inherits from the ABCMeta metaclass as _login is an abstract
    method.

    """
    def __new__(cls, name, bases, attrs): # -> Self@LoginABCMeta:
        ...
    


class BaseQuery(metaclass=LoginABCMeta):
    """
    This is the base class for all the query classes in astroquery. It
    is implemented as an abstract class and must not be directly instantiated.
    """
    def __init__(self) -> None:
        ...
    
    def __call__(self, *args, **kwargs): # -> Self@BaseQuery:
        """ init a fresh copy of self """
        ...
    


class suspend_cache:
    """
    A context manager that suspends caching.
    """
    def __init__(self, obj) -> None:
        ...
    
    def __enter__(self): # -> None:
        ...
    
    def __exit__(self, exc_type, exc_value, traceback): # -> Literal[False]:
        ...
    


class QueryWithLogin(BaseQuery):
    """
    This is the base class for all the query classes which are required to
    have a login to access the data.

    The abstract method _login() must be implemented. It is wrapped by the
    login() method, which turns off the cache. This way, login credentials
    are not stored in the cache.
    """
    def __init__(self) -> None:
        ...
    
    def login(self, *args, **kwargs): # -> None:
        ...
    
    def authenticated(self): # -> bool | None:
        ...
    


