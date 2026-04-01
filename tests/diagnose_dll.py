"""
Comprehensive ARM64 .pyd load diagnostic.
Run this AS the test script in cibuildwheel.
"""
import ctypes
import ctypes.wintypes
import os
import sys
import glob
import struct
import importlib.util
import subprocess

print("=" * 60)
print("SECTION 1: Environment")
print("=" * 60)
print(f"Python executable: {sys.executable}")
print(f"Platform tag: {__import__('sysconfig').get_platform()}")
print(f"Pointer size: {struct.calcsize('P') * 8} bit")
print(f"sys.version: {sys.version}")
print()

# Check if this Python process is actually running under emulation
print("Checking for x64 emulation...")
kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
try:
    IsWow64Process2 = kernel32.IsWow64Process2
    processMachine = ctypes.wintypes.USHORT()
    nativeMachine = ctypes.wintypes.USHORT()
    handle = kernel32.GetCurrentProcess()
    result = IsWow64Process2(handle, ctypes.byref(processMachine), ctypes.byref(nativeMachine))
    machine_names = {0: "NATIVE", 0x8664: "x86_64", 0xAA64: "ARM64", 0x14C: "x86", 0x1C4: "ARM"}
    print(f"  Process machine: {hex(processMachine.value)} ({machine_names.get(processMachine.value, 'unknown')})")
    print(f"  Native machine:  {hex(nativeMachine.value)} ({machine_names.get(nativeMachine.value, 'unknown')})")
    if processMachine.value == 0:
        print("  -> Process is running NATIVELY (no emulation)")
    else:
        print(f"  -> Process is running UNDER EMULATION!")
except Exception as e:
    print(f"  IsWow64Process2 not available: {e}")
print()

print("=" * 60)
print("SECTION 2: Test loading a known-good ARM64 extension")
print("=" * 60)
# Try importing a simple C extension that comes with Python
try:
    import _struct
    print(f"_struct loaded OK from: {getattr(_struct, '__file__', 'builtin')}")
except Exception as e:
    print(f"_struct FAILED: {e}")

try:
    import _ctypes
    print(f"_ctypes loaded OK from: {_ctypes.__file__}")
except Exception as e:
    print(f"_ctypes FAILED: {e}")

# Try a pip-installed native extension if available
for modname in ["psutil", "numpy"]:
    try:
        mod = __import__(modname)
        so_file = getattr(mod, "__file__", None)
        print(f"{modname} loaded OK from: {so_file}")
    except ImportError:
        print(f"{modname} not installed, skipping")
    except Exception as e:
        print(f"{modname} FAILED: {e}")
print()

print("=" * 60)
print("SECTION 3: Inspect the optvl .pyd binary")
print("=" * 60)
spec = importlib.util.find_spec("optvl")
if not spec:
    print("ERROR: optvl not installed!")
    sys.exit(1)

pkg_dir = spec.submodule_search_locations[0]
libs_dir = os.path.abspath(os.path.join(pkg_dir, os.pardir, "optvl.libs"))
pyd_files = glob.glob(os.path.join(pkg_dir, "libavl*.pyd"))
pyd = pyd_files[0] if pyd_files else None

print(f"Package dir: {pkg_dir}")
print(f"PYD file:    {pyd}")
print(f"PYD size:    {os.path.getsize(pyd)} bytes")
print()

# Read PE header manually for extra detail
with open(pyd, "rb") as f:
    # DOS header
    dos_magic = f.read(2)
    print(f"DOS magic: {dos_magic}")
    f.seek(0x3C)
    pe_offset = struct.unpack("<I", f.read(4))[0]
    print(f"PE offset: {pe_offset}")

    # PE signature
    f.seek(pe_offset)
    pe_sig = f.read(4)
    print(f"PE signature: {pe_sig}")

    # COFF header
    machine = struct.unpack("<H", f.read(2))[0]
    num_sections = struct.unpack("<H", f.read(2))[0]
    timestamp = struct.unpack("<I", f.read(4))[0]
    machine_names = {0x8664: "AMD64", 0xAA64: "ARM64", 0x14C: "i386"}
    print(f"Machine: {hex(machine)} ({machine_names.get(machine, 'UNKNOWN')})")
    print(f"Sections: {num_sections}")
    print(f"Timestamp: {timestamp}")

    # Optional header magic
    f.seek(pe_offset + 24)
    opt_magic = struct.unpack("<H", f.read(2))[0]
    print(f"Optional header magic: {hex(opt_magic)} ({'PE32+' if opt_magic == 0x20b else 'PE32' if opt_magic == 0x10b else 'unknown'})")
print()

print("=" * 60)
print("SECTION 4: Check Fortran runtime / hidden dependencies")
print("=" * 60)
# Use dumpbin if available, or just search for DLL strings in the binary
try:
    result = subprocess.run(
        ["dumpbin", "/dependents", pyd],
        capture_output=True, text=True, timeout=10
    )
    if result.returncode == 0:
        print("dumpbin /dependents output:")
        print(result.stdout)
    else:
        print("dumpbin failed, trying alternative...")
        raise FileNotFoundError
except (FileNotFoundError, subprocess.TimeoutExpired):
    print("dumpbin not available, scanning binary for DLL name strings...")
    with open(pyd, "rb") as f:
        data = f.read()
    # Look for common DLL name patterns
    import re
    dll_names = set()
    for match in re.finditer(rb'[\x20-\x7e]{4,}\.dll', data, re.IGNORECASE):
        name = match.group().decode('ascii', errors='ignore')
        # Filter to likely real DLL names
        if not any(c in name for c in ' <>"') and len(name) < 100:
            dll_names.add(name)
    print("DLL-like strings found in binary:")
    for name in sorted(dll_names):
        print(f"  {name}")
print()

print("=" * 60)
print("SECTION 5: DLL load tests")
print("=" * 60)

# 5a: Try loading vendored DLL directly
vendored = os.path.join(libs_dir, "msvcp140-5f1c5dd31916990d94181e07bc3afb32.dll")
print(f"Libs dir exists: {os.path.isdir(libs_dir)}")
if os.path.isdir(libs_dir):
    print(f"Libs contents: {os.listdir(libs_dir)}")
print(f"Vendored DLL exists: {os.path.exists(vendored)}")

if os.path.exists(vendored):
    try:
        h = ctypes.WinDLL(vendored)
        print(f"5a - Direct load vendored msvcp140: SUCCESS")
    except OSError as e:
        print(f"5a - Direct load vendored msvcp140: FAILED ({e})")
        # If even the vendored DLL can't load, check ITS deps
        print("     This DLL depends on VCRUNTIME140.dll - checking that...")
print()

# 5b: Load .pyd WITHOUT os.add_dll_directory
print("5b - Load .pyd WITHOUT os.add_dll_directory:")
try:
    h = ctypes.WinDLL(pyd)
    print("  SUCCESS")
except OSError as e:
    print(f"  FAILED: {e}")
print()

# 5c: Add dll directory and retry
print(f"5c - Adding os.add_dll_directory('{libs_dir}')")
try:
    os.add_dll_directory(libs_dir)
    print("  os.add_dll_directory succeeded")
except OSError as e:
    print(f"  os.add_dll_directory FAILED: {e}")

try:
    h = ctypes.WinDLL(pyd)
    print("  Load .pyd AFTER add_dll_directory: SUCCESS")
except OSError as e:
    print(f"  Load .pyd AFTER add_dll_directory: FAILED ({e})")
print()

# 5d: Copy vendored DLL next to .pyd and try
print("5d - Copy vendored DLL next to .pyd and try:")
import shutil
import tempfile
tmpdir = tempfile.mkdtemp()
pyd_copy = shutil.copy(pyd, tmpdir)
if os.path.exists(vendored):
    shutil.copy(vendored, tmpdir)
    # Also rename to unmangled name
    shutil.copy(vendored, os.path.join(tmpdir, "msvcp140.dll"))
print(f"  Temp dir contents: {os.listdir(tmpdir)}")
try:
    h = ctypes.WinDLL(pyd_copy)
    print("  Load .pyd with vendored DLL alongside: SUCCESS")
except OSError as e:
    print(f"  Load .pyd with vendored DLL alongside: FAILED ({e})")
print()

# 5e: Try LoadLibraryEx with LOAD_WITH_ALTERED_SEARCH_PATH
print("5e - LoadLibraryExW with LOAD_WITH_ALTERED_SEARCH_PATH:")
LoadLibraryExW = kernel32.LoadLibraryExW
LoadLibraryExW.argtypes = [ctypes.wintypes.LPCWSTR, ctypes.wintypes.HANDLE, ctypes.wintypes.DWORD]
LoadLibraryExW.restype = ctypes.wintypes.HMODULE
LOAD_WITH_ALTERED_SEARCH_PATH = 0x00000008
h = LoadLibraryExW(pyd, None, LOAD_WITH_ALTERED_SEARCH_PATH)
if h:
    print(f"  SUCCESS (handle={h})")
else:
    err = ctypes.get_last_error()
    print(f"  FAILED: Win32 error {err}")
print()

# 5f: Try with LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR | LOAD_LIBRARY_SEARCH_DEFAULT_DIRS
print("5f - LoadLibraryExW with SEARCH_DLL_LOAD_DIR | SEARCH_DEFAULT_DIRS:")
LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR = 0x00000100
LOAD_LIBRARY_SEARCH_DEFAULT_DIRS = 0x00001000
h = LoadLibraryExW(pyd, None, LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR | LOAD_LIBRARY_SEARCH_DEFAULT_DIRS)
if h:
    print(f"  SUCCESS (handle={h})")
else:
    err = ctypes.get_last_error()
    print(f"  FAILED: Win32 error {err}")
print()

print("=" * 60)
print("SECTION 6: Full import test")
print("=" * 60)
try:
    import optvl
    print("import optvl: SUCCESS")
    print(f"optvl.__file__: {optvl.__file__}")
except Exception as e:
    print(f"import optvl: FAILED ({e})")

try:
    from optvl import OVLSolver
    print("from optvl import OVLSolver: SUCCESS")
except Exception as e:
    print(f"from optvl import OVLSolver: FAILED ({e})")
print()
print("Done.")