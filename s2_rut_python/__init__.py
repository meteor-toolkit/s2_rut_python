"""s2_rut_python - Pure Python port of SNAP's Sentinel-2 radiometric uncertainty tool"""

from importlib.metadata import version as _pkg_version, PackageNotFoundError

__author__ = "MetEOR Toolkit Team"
__all__: list = []

try:
    __version__ = _pkg_version("s2_rut_python")
except PackageNotFoundError:
    __version__ = "0.0.0"
