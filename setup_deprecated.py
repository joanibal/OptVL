from setuptools import setup
import re
import os
import subprocess
import re
from pathlib import Path


this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()


def make_avl_lib():
    cwd = os.getcwd()
    # os.chdir("optvl")
    # subprocess.run(["make", "clean"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    config_file = "config/defaults/config.LINUX_GFORTRAN.mk"

    build_dir = os.path.join(cwd, "src/build")
    subprocess.run(["cp", config_file, "config/config.mk"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    # subprocess.run(["make", "clean"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    p2 = subprocess.run(["make"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    compile_log = os.path.join(build_dir, "compile.log")
    with open(compile_log, "wb") as f:
        f.write(p2.stdout)
    if p2.returncode != 0:
        with open(compile_log, "r") as f:
            print(f.read())
        raise OSError("make", f"The compile command failed! Check the log at {compile_log} for more information.")


def get_version_from_pyproject():
    pyproject = Path(__file__).parent / "pyproject.toml"
    text = pyproject.read_text()

    match = re.search(r'^version\s*=\s*"([^"]+)"', text, re.MULTILINE)
    if match:
        return match.group(1)
    raise RuntimeError("Version not found in pyproject.toml")


if __name__ == "__main__":
    # to install locally use...
    # `python setup_deprecated.py develop`

    print(get_version_from_pyproject())

    setup(
        name="optvl",
        version=get_version_from_pyproject(),
        description="A direct Python interface for Mark Drela and Harold Youngren's famous AVL code.",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://github.com/joanibal/optvl",
        license="GPL-2.0",
        packages=[
            "optvl",
        ],
        package_data={"optvl": ["*.so"]},
        install_requires=["numpy"],
        extras_require={"plotting": ["matplotlib"], "testing": ["testflo>=1.4.5"]},
        classifiers=["Programming Language :: Python, Fortran"],
    )
