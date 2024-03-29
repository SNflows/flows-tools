"""
Helper script for getting the brightest star in a Hawk-I (or other)
instrument Field of View for a given FLOWS target, retrieved via
tendrils API. Extensible for other instruments and projects.
Can optionally make a finding chart using NASA SkyView.
"""
from . import make_fc as make_fc  # noqa: F401
from . import run_get_brightest as get_brightest  # noqa: F401

# from .run_get_brightest import main as main  # noqa: F401
from .version import __version__ as __version__  # noqa: F401
