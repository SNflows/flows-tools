"""
This type stub file was generated by pyright.
"""

import configparser
import requests
from dataclasses import dataclass
from typing import Optional

"""
For setting up and accessing remote urls. Defaults to FLOWS project ones but you can set up your own.
"""

@dataclass
class RemoteUrls:
    base_url: str = ...
    datafiles_url: str = ...
    targets_url: str = ...
    sites_url: str = ...
    set_photometry_status_url: str = ...
    photometry_url: str = ...
    photometry_upload_url: str = ...
    cleanup_photometry_status_url: str = ...
    catalogs_url: str = ...
    catalogs_missing_url: str = ...
    filters_url: str = ...
    lightcurves_url: str = ...
    targets_post_url: str = ...
    @staticmethod
    def urls_from_config(config: Optional[configparser.ConfigParser] = ...) -> dict:
        """
        Update default remote urls with values found from config. Will warn if ALL values are not found.
        Returns: dict[str:str] representing urls
        These should be under a [URLS] section in config.ini or in another config.
        """
        ...
    def update(self, new: Optional[dict] = ...) -> None: ...
    def make_convenient(self) -> None: ...

def urls(urls_instance: Optional[RemoteUrls] = ...) -> RemoteUrls:
    """
    Convenience function for getting an RemoteUrls instance.
    Warning! Destructively overrides all urls that are not "base_url"
    Returns: RemoteUrls , but with base_url prepended.
    """
    ...

def urls_from_config(config: Optional[configparser.ConfigParser] = ...) -> RemoteUrls:
    """
    Use config to return RemoteUrls instance with values read from default config or given ConfigParser instance.
    returned RemoteUrls have base_url prepended, similar to urls().
    Args:
        config: Optional[ConfigParser] instance which defaults to None. Uses default config if None.

    Returns: RemoteUrls instance with the urls populated from the given config
    """
    ...

def get_request(
    url: str, token: str = ..., params: Optional[dict] = ..., headers: Optional[dict] = ...
) -> requests.Response:
    """
    Make a get request using request with given url, params, and header. Token can be given inplace of header
    to create a default header that uses the token.
    Args:
        url: str = url for get request
        token: Optional[str] = api token. Can also be provided as part of headers dict.
        params: Optional[dict] = dict of params
        headers: Optional[dict] = headers dict. Can also be created from just the token.

    Returns: requests.Response

    """
    ...

def post_request(
    url: str,
    token: str = ...,
    params: Optional[dict] = ...,
    data: Optional[dict] = ...,
    headers: Optional[dict] = ...,
    files: Optional[dict] = ...,
) -> requests.Response:
    """
    Make a post request using request with given url, params, and header. Token can be given inplace of header
    to create a default header that uses the token.
    Args:
        url: str = url for get request
        token: str, optional = api token. Can also be provided as part of headers dict.
        params: dict, optional = dict of params
        data: dict, optional = dict of data
        headers: dict, optional = headers dict. Can also be created from just the token.
        files: dict, optional = files dict.

    Returns: requests.Response
    """
    ...
