"""
Helper script for getting the brightest star in a Hawk-I (or other)
instrument Field of View for a given FLOWS target, retrieved via
tendrils API. Extensible for other instruments and projects.
Can optionally make a finding chart using NASA SkyView.
"""
from .version import __version__
from .run_get_brightest import main