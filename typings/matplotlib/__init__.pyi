from _typeshed import Incomplete
from collections.abc import Generator, MutableMapping
from matplotlib._api import MatplotlibDeprecationWarning as MatplotlibDeprecationWarning
from matplotlib.cbook import sanitize_sequence as sanitize_sequence
from matplotlib.rcsetup import cycler as cycler, validate_backend as validate_backend
from typing import NamedTuple

__bibtex__: str

class _VersionInfo(NamedTuple):
    major: Incomplete
    minor: Incomplete
    micro: Incomplete
    releaselevel: Incomplete
    serial: Incomplete

class __getattr__:
    __version_info__: Incomplete
    URL_REGEX: Incomplete

def set_loglevel(level) -> None: ...

class _ExecInfo(NamedTuple):
    executable: Incomplete
    raw_version: Incomplete
    version: Incomplete

class ExecutableNotFoundError(FileNotFoundError): ...

def checkdep_usetex(s): ...
def get_configdir(): ...
def get_cachedir(): ...
def get_data_path(): ...
def matplotlib_fname(): ...

class RcParams(dict, MutableMapping):
    validate: Incomplete
    def __init__(self, *args, **kwargs) -> None: ...
    def __setitem__(self, key, val) -> None: ...
    def __getitem__(self, key): ...
    def __iter__(self): ...
    def __len__(self) -> int: ...
    def find_all(self, pattern): ...
    def copy(self): ...

def rc_params(fail_on_error: bool = ...): ...
def is_url(filename): ...
def rc_params_from_file(fname, fail_on_error: bool = ..., use_default_template: bool = ...): ...

rcParamsDefault: Incomplete
rcParams: Incomplete
rcParamsOrig: Incomplete
defaultParams: Incomplete

def rc(group, **kwargs) -> None: ...
def rcdefaults() -> None: ...
def rc_file_defaults() -> None: ...
def rc_file(fname, *, use_default_template: bool = ...) -> None: ...
def rc_context(rc: Incomplete | None = ..., fname: Incomplete | None = ...) -> Generator[None, None, None]: ...
def use(backend, *, force: bool = ...) -> None: ...
def get_backend(): ...
def interactive(b) -> None: ...
def is_interactive(): ...

default_test_modules: Incomplete

def test(verbosity: Incomplete | None = ..., coverage: bool = ..., **kwargs): ...
