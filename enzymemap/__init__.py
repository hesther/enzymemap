"""Python package to atom-map, correct and suggest enzymatic reactions"""

# Add imports here
from .enzymemap import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
