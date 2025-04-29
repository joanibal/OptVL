import importlib.metadata

try:
    # __package__ allows for the case where __name__ is "__main__"
    __version__ = importlib.metadata.version(__package__ or __name__)
except importlib.metadata.PackageNotFoundError:
    __version__ = "0.0.0"

from .optvl_class import OVLSolver

try:
    from .om_wrapper import AVLGroup, AVLMeshReader, Differencer
except ImportError:
    # if openmdao is not installed, then we can't use the wrapper
    pass
