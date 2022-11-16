"""
This type stub file was generated by pyright.
"""

from astropy.utils.exceptions import AstropyWarning

"""
Custom exceptions used in the astroquery query classes
"""
__all__ = ['TimeoutError', 'InvalidQueryError', 'RemoteServiceError', 'TableParseError', 'LoginError', 'ResolverError', 'NoResultsWarning', 'LargeQueryWarning', 'InputWarning', 'AuthenticationWarning', 'MaxResultsWarning', 'CorruptDataWarning']
class TimeoutError(Exception):
    """
    Raised on failure to establish connection with server
    within a particular time limit
    """
    ...


class InvalidQueryError(Exception):
    """
    Errors related to invalid queries.
    """
    ...


class TableParseError(Exception):
    """
    Errors related to VOTable parsing.
    These should be either submitted as issues to astropy or to the originating
    service.
    """
    ...


class RemoteServiceError(Exception):
    """
    Errors related to the remote service, i.e. if the service returns an error
    page
    """
    ...


class LoginError(Exception):
    """
    Errors due to failed logins.  Should only be raised for services for which
    a login is a prerequisite for the requested action
    """
    ...


class ResolverError(Exception):
    """
    Errors due to failing to resolve an object name/id to a specific
    sky coordinate.
    """
    ...


class NoResultsWarning(AstropyWarning):
    """
    Astroquery warning class to be issued when a query returns no result.
    """
    ...


class LargeQueryWarning(AstropyWarning):
    """
    Astroquery warning class to be issued when a query is larger than
    recommended for a given service.
    """
    ...


class InputWarning(AstropyWarning):
    """
    Astroquery warning class to be issued when user input is incorrect
    in some way but doesn't prevent the function from running.
    """
    ...


class AuthenticationWarning(AstropyWarning):
    """
    Astroquery warning class to be issued when there are problems with
    user authentication.
    """
    ...


class MaxResultsWarning(AstropyWarning):
    """
    Astroquery warning class to be issued when the maximum allowed
    results are returned.
    """
    ...


class CorruptDataWarning(AstropyWarning):
    """
    Astroquery warning class to be issued when there is a sign that the
    (partially) downloaded data are corrupt.
    """
    ...


class EmptyResponseError(ValueError):
    """
    Astroquery error class to be raised when the query returns an empty result
    """
    ...


