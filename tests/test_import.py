# =============================================================================
# Extension modules
# =============================================================================
import struct, sys, sysconfig
print("Pointer size:", struct.calcsize("P") * 8, "bit")  # 64 for both, less useful
print("Platform tag:", sysconfig.get_platform())          # THIS is the key one
print("sys.executable:", sys.executable)

import ctypes
import ctypes.wintypes
import os, sys, struct, glob

print("Python arch:", struct.calcsize("P") * 8, "bit")
print("Executable:", sys.executable)

# Check what VCRUNTIME140.dll Python itself loaded
import ctypes.util
kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
GetModuleFileNameW = kernel32.GetModuleFileNameW
GetModuleFileNameW.argtypes = [ctypes.wintypes.HMODULE, ctypes.wintypes.LPWSTR, ctypes.wintypes.DWORD]
GetModuleFileNameW.restype = ctypes.wintypes.DWORD
GetModuleHandleW = kernel32.GetModuleHandleW
GetModuleHandleW.argtypes = [ctypes.wintypes.LPCWSTR]
GetModuleHandleW.restype = ctypes.wintypes.HMODULE

for dll_name in ["VCRUNTIME140.dll", "msvcp140.dll", "python312.dll"]:
    h = GetModuleHandleW(dll_name)
    if h:
        buf = ctypes.create_unicode_buffer(260)
        GetModuleFileNameW(h, buf, 260)
        print(f"{dll_name}: loaded from {buf.value}")
    else:
        print(f"{dll_name}: not loaded")

# Search common locations for VCRUNTIME140.dll
search_dirs = [
    os.path.dirname(sys.executable),
    os.path.join(os.environ.get("SystemRoot", r"C:\Windows"), "System32"),
    os.path.join(os.environ.get("SystemRoot", r"C:\Windows"), "SysWOW64"),
]
for d in search_dirs:
    path = os.path.join(d, "VCRUNTIME140.dll")
    if os.path.exists(path):
        # Read PE machine type
        with open(path, "rb") as f:
            f.seek(0x3C)  # e_lfanew offset
            pe_offset = int.from_bytes(f.read(4), "little")
            f.seek(pe_offset + 4)  # Machine field
            machine = int.from_bytes(f.read(2), "little")
        arch = {0x8664: "x86_64", 0xAA64: "ARM64", 0x14C: "x86"}.get(machine, hex(machine))
        print(f"Found {path} -> {arch}")

# Try a direct load of the .pyd to get better error info
print("\n--- Direct load test ---")
import importlib.util
spec = importlib.util.find_spec("optvl")
if spec:
    pyd = glob.glob(os.path.join(spec.submodule_search_locations[0], "libavl*.pyd"))
    if pyd:
        print(f"Attempting ctypes load of: {pyd[0]}")
        try:
            ctypes.WinDLL(pyd[0])
            print("SUCCESS: Direct load works")
        except OSError as e:
            print(f"FAILED: {e}")
        
        # Now test loading from a temp copy (what MExt does)
        import shutil, tempfile
        tmpdir = tempfile.mkdtemp()
        copied = shutil.copy(pyd[0], tmpdir)
        print(f"Attempting ctypes load of temp copy: {copied}")
        try:
            ctypes.WinDLL(copied)
            print("SUCCESS: Temp copy load works")
        except OSError as e:
            print(f"FAILED from temp: {e}")
# =============================================================================
# Standard Python Modules
# =============================================================================
import os

# =============================================================================
# External Python modules
# =============================================================================
import unittest

base_dir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder
geom_dir = os.path.join(base_dir, '..', 'geom_files')

geom_file = os.path.join(geom_dir, "aircraft.avl")
mass_file = os.path.join(geom_dir, "aircraft.mass")


class TestImport(unittest.TestCase):
    # TODO: add test for expected input output errors
    def test_instances(self):
        from optvl import OVLSolver

        ovl_solver1 = OVLSolver(geo_file=geom_file)

        ovl_solver2 = OVLSolver(geo_file=geom_file)

        assert ovl_solver1.avl is not ovl_solver2.avl

        ovl_solver1.set_avl_fort_arr("CASE_R", "ALFA", 1.1)
        ovl_solver2.set_avl_fort_arr("CASE_R", "ALFA", 2.0)

        assert ovl_solver1.get_avl_fort_arr("CASE_R", "ALFA") != ovl_solver2.get_avl_fort_arr("CASE_R", "ALFA")


if __name__ == "__main__":
    unittest.main()
