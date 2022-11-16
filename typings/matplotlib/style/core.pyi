from _typeshed import Incomplete
from collections.abc import Generator

class __getattr__:
    STYLE_FILE_PATTERN: Incomplete

def use(style): ...
def context(style, after_reset: bool = ...) -> Generator[None, None, None]: ...

class _StyleLibrary(dict):
    def __getitem__(self, key): ...

library: Incomplete
available: Incomplete

def reload_library() -> None: ...
