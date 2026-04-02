"""s2_rut_python - Pure Python port of SNAP's Sentinel-2 radiometric uncertainty tool"""

__author__ = "Sam Hunt <sam.hunt@npl.co.uk>"
__all__: list = []

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions
