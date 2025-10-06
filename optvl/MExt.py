"""
This file was pulled from ADflow on March 11, 2022.
"""

import tempfile
from importlib.util import find_spec
from pathlib import Path
import os
import shutil
import sys
import platform


def _tmp_pkg(tempDir):
    """
    Create a temporary package.

    Returns (name, path)
    """
    while True:
        path = tempfile.mkdtemp(dir=tempDir)
        name = os.path.basename(path)
        spec = find_spec(name)
        # None means the name was not found
        if spec is None:
            break
        # if name is found, delete and try again
        os.rmdir(path)
    # this creates an init file so that python recognizes this as a package
    Path(os.path.join(path, "__init__.py")).touch()
    return name, path


class MExt(object):
    """
    Load a unique copy of a module that can be treated as a "class instance".
    """

    def __init__(self, libName, packageName, pip_name, lib_so_file=None, debug=False):
        if lib_so_file is None:
            lib_so_file = f"{libName}.so"

        tmpdir = tempfile.gettempdir()
        self.name = libName
        self.debug = debug
        # first find the "real" module on the "real" syspath
        spec = find_spec(packageName)
        srcpath = os.path.join(spec.submodule_search_locations[0], lib_so_file)
        # now create a temp directory for the bogus package
        self._pkgname, self._pkgdir = _tmp_pkg(tmpdir)
        # copy the original module to the new package
        shutil.copy(srcpath, self._pkgdir)
        if platform.system() == "Darwin":
            # create a sym link to the orginal module .dylibs folder
            blas_libs_dir = ".dylibs"
            source_path = os.path.join(spec.submodule_search_locations[0], blas_libs_dir)
            target_path = os.path.join(self._pkgdir, blas_libs_dir)

            if not os.path.exists(target_path) and os.path.exists(source_path):
                # print("Creating symlink from {} to {}".format(source_path, target_path))
                os.symlink(source_path, target_path)

        elif platform.system() == "Linux":
            blas_libs_dir = f"{pip_name}.libs"
            source_path = os.path.join(spec.submodule_search_locations[0], "..", blas_libs_dir)
            target_path = os.path.join("/tmp", blas_libs_dir)

            if not os.path.exists(target_path) and os.path.exists(source_path):
                # print("Creating symlink from {} to {}".format(source_path, target_path))
                os.symlink(source_path, target_path)
        else:
            # raise NotImplementedError("platform not supported")
            print(tmpdir)
            srcpath = os.path.join(spec.submodule_search_locations[0], "libavl.cp39-win_amd64.dll.a")
            pass
        # add the directory containing the new package to the search path
        sys.path.append(tmpdir)
        # import the module
        # __import__ returns the package, not the sub-module

        self._pkg = __import__(self._pkgname, globals(), locals(), [self.name])
        # remove the bogus directory from sys.path
        sys.path.remove(tmpdir)
        # return the module object
        self._module = getattr(self._pkg, self.name)
        # now add the module's stuff to this class
        self.__dict__.update(self._module.__dict__)

    def __del__(self):
        # remove module if not in debug mode
        if not self.debug:
            # if the module was imported, remove it from sys.modules
            if hasattr(self, "_pkg"):
                del sys.modules[self._module.__name__]
                del sys.modules[self._pkg.__name__]

            # now try to delete the files and directory
            # shutil.rmtree(self._pkgdir)
            # make sure the original module is loaded -
            # otherwise python crashes on exit
            # if MExt objects have not been explicitly 'del'd,
            # and __del__ is occurring on python shutdown, the import will fail
            # and the exception is caught here
            try:
                __import__(self.name)
            except ImportError:
                pass
