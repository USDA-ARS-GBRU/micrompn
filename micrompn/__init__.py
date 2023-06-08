from importlib.metadata import version, PackageNotFoundError

from .mpn import *
from .main import *



__all__ = ["mpn", "main"]


try:
    __version__ = version("micrompn")
except PackageNotFoundError:
    # package is not installed
    pass