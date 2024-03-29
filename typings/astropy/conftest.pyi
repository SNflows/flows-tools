"""
This type stub file was generated by pyright.
"""

import sys
import pytest
from astropy.utils.compat.optional_deps import HAS_MATPLOTLIB

"""
This file contains pytest configuration settings that are astropy-specific
(i.e.  those that would not necessarily be shared by affiliated packages
making use of astropy's test runner).
"""
if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
    ...
if HAS_MATPLOTLIB:
    ...
matplotlibrc_cache = ...
@pytest.fixture
def ignore_matplotlibrc(): # -> Generator[None, None, None]:
    ...

@pytest.fixture
def fast_thread_switching(): # -> Generator[None, None, None]:
    """Fixture that reduces thread switching interval.

    This makes it easier to provoke race conditions.
    """
    ...

def pytest_configure(config): # -> None:
    ...

def pytest_unconfigure(config): # -> None:
    ...

def pytest_terminal_summary(terminalreporter): # -> None:
    """Output a warning to IPython users in case any tests failed."""
    ...

