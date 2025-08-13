from setuptools import setup
import re
import os


this_directory = os.path.abspath(os.path.dirname(__file__))
project_directory = os.path.dirname(this_directory)
optvl_directory  = os.path.join(project_directory, 'optvl')
sym_link = os.path.join(this_directory, 'optvl')

with open(os.path.join(project_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()
    

if not os.path.exists(sym_link):
    # Try to create a relative symlink (POSIX-friendly)
    try:
        os.symlink(optvl_directory, sym_link)
    except (OSError, NotImplementedError) as exc:
        # On Windows or if not permitted, fall through to wrapper creation
        # print a helpful message only once
        print(f"Could not create symlink {sym_link} -> {optvl_directory}: {exc}\n")
        exit
        

def get_version_from_pyproject():
    path  = os.path.join(project_directory, "pyproject.toml")
    with open(path, 'r', encoding='utf-8') as f:
        text = f.read()

    match = re.search(r'^version\s*=\s*"([^"]+)"', text, re.MULTILINE)
    if match:
        return match.group(1)
    raise RuntimeError("Version not found in pyproject.toml")


if __name__ == "__main__":
    # after running make in the top
    # `pip install`


    setup(
        name="optvl",
        version=f"{get_version_from_pyproject()}.dev0",
        description="A direct Python interface for Mark Drela and Harold Youngren's famous AVL code.",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://github.com/joanibal/optvl",
        license="GPL-2.0",
        packages=["optvl"],
        package_dir={"optvl": "../"},
        package_data={"optvl": ["*.so"]},
        install_requires=["numpy"],
        extras_require={"plotting": ["matplotlib"], "testing": ["testflo>=1.4.5"]},
        classifiers=["Programming Language :: Python, Fortran"],
    )
