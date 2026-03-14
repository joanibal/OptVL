# from scipy

set -xe

WHEEL="$1"
DEST_DIR="$2"


# Skip the strip command based on TARGET_ARCH
# TARGET_ARCH should be set by the CI environment (e.g., ARM64, AMD64)
TARGET_ARCH="${TARGET_ARCH:-}" # Default to empty string if not set

if [ "$TARGET_ARCH" = "ARM64" ]; then
    echo "Skipping stripping for ARM64 target."
    
    pushd $DEST_DIR
    mkdir -p tmp
    pushd tmp
    wheel unpack $WHEEL
    pushd optvl*
    ls 
    llvm-objdump -p *.pyd | grep 'DLL Name'
    dumpbin /dependents *.pyd
else
    echo "Performing stripping for AMD64 target."

    # create a temporary directory in the destination folder and unpack the wheel
    # into there
    pushd $DEST_DIR
    mkdir -p tmp
    pushd tmp
    wheel unpack $WHEEL
    pushd optvl*

    # To avoid DLL hell, the file name of libopenblas that's being vendored with
    # the wheel has to be name-mangled. delvewheel is unable to name-mangle PYD
    # containing extra data at the end of the binary, which frequently occurs when
    # building with mingw.
    # We therefore find each PYD in the directory structure and strip them.

    for f in $(find ./optvl* -name '*.pyd'); do strip $f; done


    # now repack the wheel and overwrite the original
    wheel pack .
    mv -fv *.whl $WHEEL

    cd $DEST_DIR
    rm -rf tmp
fi

delvewheel repair  -w $DEST_DIR $WHEEL

