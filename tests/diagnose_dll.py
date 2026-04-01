"""
Final diagnostic: What exactly is wrong with this .pyd?
Run as the test script in cibuildwheel.
"""
import ctypes
import ctypes.wintypes
import os
import sys
import glob
import struct
import importlib
import importlib.util

def read_pe_details(filepath):
    """Read detailed PE structure."""
    info = {}
    with open(filepath, "rb") as f:
        data = f.read()

    # DOS header
    pe_offset = struct.unpack_from("<I", data, 0x3C)[0]

    # COFF header (at pe_offset + 4, after "PE\0\0")
    coff_start = pe_offset + 4
    machine, num_sections, timestamp, sym_table_ptr, num_symbols, opt_header_size, characteristics = \
        struct.unpack_from("<HHIIIHH", data, coff_start)

    info["machine"] = machine
    info["num_sections"] = num_sections
    info["opt_header_size"] = opt_header_size
    info["characteristics"] = characteristics
    info["file_size"] = len(data)

    # Optional header
    opt_start = coff_start + 20
    opt_magic = struct.unpack_from("<H", data, opt_start)[0]
    info["opt_magic"] = opt_magic

    if opt_magic == 0x20b:  # PE32+
        info["image_base"] = struct.unpack_from("<Q", data, opt_start + 24)[0]
        info["section_alignment"] = struct.unpack_from("<I", data, opt_start + 32)[0]
        info["file_alignment"] = struct.unpack_from("<I", data, opt_start + 36)[0]
        info["size_of_image"] = struct.unpack_from("<I", data, opt_start + 56)[0]
        info["size_of_headers"] = struct.unpack_from("<I", data, opt_start + 60)[0]
        info["checksum"] = struct.unpack_from("<I", data, opt_start + 64)[0]
        info["subsystem"] = struct.unpack_from("<H", data, opt_start + 68)[0]
        info["dll_characteristics"] = struct.unpack_from("<H", data, opt_start + 70)[0]
        info["num_data_dirs"] = struct.unpack_from("<I", data, opt_start + 108)[0]

        # Data directories start at opt_start + 112
        dd_start = opt_start + 112
        dd_names = ["Export", "Import", "Resource", "Exception", "Certificate",
                     "BaseReloc", "Debug", "Architecture", "GlobalPtr", "TLS",
                     "LoadConfig", "BoundImport", "IAT", "DelayImport",
                     "CLRRuntime", "Reserved"]
        info["data_dirs"] = {}
        for i in range(min(info["num_data_dirs"], 16)):
            rva, size = struct.unpack_from("<II", data, dd_start + i * 8)
            if size > 0:
                info["data_dirs"][dd_names[i] if i < len(dd_names) else f"Dir{i}"] = (rva, size)

    # Section headers
    section_start = opt_start + opt_header_size
    info["sections"] = []
    for i in range(num_sections):
        s_off = section_start + i * 40
        name = data[s_off:s_off+8].rstrip(b'\x00').decode('ascii', errors='replace')
        vsize, vrva, raw_size, raw_ptr, _, _, _, _, characteristics = \
            struct.unpack_from("<IIIIIIHHI", data, s_off + 8)
        info["sections"].append({
            "name": name, "virtual_size": vsize, "virtual_rva": vrva,
            "raw_size": raw_size, "raw_ptr": raw_ptr,
            "characteristics": characteristics
        })

    # Check export table for PyInit_
    if "Export" in info.get("data_dirs", {}):
        export_rva, export_size = info["data_dirs"]["Export"]
        # Find which section contains this RVA
        for sec in info["sections"]:
            if sec["virtual_rva"] <= export_rva < sec["virtual_rva"] + sec["virtual_size"]:
                file_offset = sec["raw_ptr"] + (export_rva - sec["virtual_rva"])
                # Export directory table
                num_funcs = struct.unpack_from("<I", data, file_offset + 20)[0]
                num_names = struct.unpack_from("<I", data, file_offset + 24)[0]
                names_rva = struct.unpack_from("<I", data, file_offset + 32)[0]
                names_file_offset = sec["raw_ptr"] + (names_rva - sec["virtual_rva"])
                exports = []
                for j in range(min(num_names, 20)):
                    name_rva = struct.unpack_from("<I", data, names_file_offset + j * 4)[0]
                    name_foff = sec["raw_ptr"] + (name_rva - sec["virtual_rva"])
                    end = data.index(b'\x00', name_foff)
                    exp_name = data[name_foff:end].decode('ascii', errors='replace')
                    exports.append(exp_name)
                info["exports"] = exports
                break

    # Check LoadConfig for ARM64X / CHPE metadata
    if "LoadConfig" in info.get("data_dirs", {}):
        lc_rva, lc_size = info["data_dirs"]["LoadConfig"]
        for sec in info["sections"]:
            if sec["virtual_rva"] <= lc_rva < sec["virtual_rva"] + sec["virtual_size"]:
                lc_off = sec["raw_ptr"] + (lc_rva - sec["virtual_rva"])
                lc_struct_size = struct.unpack_from("<I", data, lc_off)[0]
                info["load_config_size"] = lc_struct_size
                # On ARM64, offset 200 (0xC8) contains CHPEMetadataPointer
                # if this exists and is nonzero, it's ARM64X/CHPE
                if lc_struct_size >= 0xCC:
                    chpe_ptr = struct.unpack_from("<Q", data, lc_off + 0xC8)[0]
                    info["chpe_metadata_pointer"] = chpe_ptr
                break

    # Quick check: sample first 16 bytes of .text section for instruction patterns
    for sec in info["sections"]:
        if sec["name"] == ".text":
            sample_start = sec["raw_ptr"]
            sample = data[sample_start:sample_start+64]
            info["text_sample_hex"] = sample.hex()
            # ARM64 instructions are always 4 bytes
            # Common ARM64 patterns: STP (store pair) often starts functions
            # x86_64 patterns: push rbp = 0x55, mov rbp,rsp = 0x48 0x89 0xe5
            first_bytes = list(sample[:8])
            info["text_first_bytes"] = first_bytes
            break

    return info


print("=" * 60)
print("SECTION 1: Python native import of the .pyd")
print("=" * 60)

# This is the critical test: can Python's own import machinery load it?
spec = importlib.util.find_spec("optvl")
pkg_dir = spec.submodule_search_locations[0]
libs_dir = os.path.abspath(os.path.join(pkg_dir, os.pardir, "optvl.libs"))

# First ensure delvewheel patch is active
if os.path.isdir(libs_dir):
    os.add_dll_directory(libs_dir)
    print(f"Added DLL directory: {libs_dir}")

# Try Python-native import of the .pyd
print("\nAttempting: import optvl.libavl")
try:
    import optvl.libavl
    print("SUCCESS - Python's import system can load the .pyd!")
except ImportError as e:
    print(f"FAILED: {e}")

# Also try importlib directly
print("\nAttempting: importlib.import_module('optvl.libavl')")
try:
    mod = importlib.import_module("optvl.libavl")
    print(f"SUCCESS: {mod}")
except ImportError as e:
    print(f"FAILED: {e}")
print()

print("=" * 60)
print("SECTION 2: Compare PE structure - broken .pyd vs working .pyd")
print("=" * 60)

pyd_files = glob.glob(os.path.join(pkg_dir, "libavl*.pyd"))
broken_pyd = pyd_files[0] if pyd_files else None

# Find a known-working .pyd for comparison
working_pyd = None
for candidate in [
    os.path.join(os.path.dirname(ctypes.__file__), "_ctypes.pyd"),  # might not exist here
]:
    if os.path.exists(candidate):
        working_pyd = candidate
        break

# Try _ctypes from DLLs directory
if not working_pyd:
    import _ctypes
    working_pyd = _ctypes.__file__

# Also try psutil
psutil_pyd = None
try:
    import psutil._psutil_windows as pw
    psutil_pyd = pw.__file__
except:
    pass

print(f"Broken .pyd:  {broken_pyd}")
print(f"Working .pyd: {working_pyd}")
if psutil_pyd:
    print(f"Psutil .pyd:  {psutil_pyd}")
print()

targets = [("BROKEN (libavl)", broken_pyd), ("WORKING (_ctypes)", working_pyd)]
if psutil_pyd:
    targets.append(("WORKING (psutil)", psutil_pyd))

machine_names = {0x8664: "AMD64", 0xAA64: "ARM64", 0x14C: "i386", 0xA641: "ARM64EC"}

for label, path in targets:
    if not path or not os.path.exists(path):
        print(f"\n--- {label}: FILE NOT FOUND ---")
        continue
    print(f"\n--- {label}: {os.path.basename(path)} ---")
    try:
        info = read_pe_details(path)
        print(f"  Machine:            {hex(info['machine'])} ({machine_names.get(info['machine'], 'UNKNOWN')})")
        print(f"  File size:          {info['file_size']}")
        print(f"  Opt magic:          {hex(info['opt_magic'])}")
        print(f"  Characteristics:    {hex(info['characteristics'])}")
        print(f"  DLL chars:          {hex(info.get('dll_characteristics', 0))}")
        print(f"  Subsystem:          {info.get('subsystem', '?')}")
        print(f"  Section alignment:  {info.get('section_alignment', '?')}")
        print(f"  File alignment:     {info.get('file_alignment', '?')}")
        print(f"  Image size:         {info.get('size_of_image', '?')}")
        print(f"  Header size:        {info.get('size_of_headers', '?')}")
        print(f"  Num data dirs:      {info.get('num_data_dirs', '?')}")
        print(f"  Load config size:   {info.get('load_config_size', 'N/A')}")
        print(f"  CHPE metadata ptr:  {hex(info.get('chpe_metadata_pointer', 0))}")

        print(f"  Data directories:")
        for name, (rva, size) in info.get("data_dirs", {}).items():
            print(f"    {name:20s} RVA={hex(rva):10s} Size={size}")

        print(f"  Sections:")
        for sec in info["sections"]:
            print(f"    {sec['name']:8s} VSize={sec['virtual_size']:8d} VRVA={hex(sec['virtual_rva']):10s} "
                  f"RawSize={sec['raw_size']:8d} RawPtr={hex(sec['raw_ptr']):10s} "
                  f"Chars={hex(sec['characteristics'])}")

        if "exports" in info:
            print(f"  Exports (first 20):")
            for exp in info["exports"]:
                print(f"    {exp}")

        if "text_first_bytes" in info:
            print(f"  .text first 8 bytes: {info['text_first_bytes']}")
            print(f"  .text first 64B hex: {info.get('text_sample_hex', 'N/A')}")

    except Exception as e:
        print(f"  ERROR reading PE: {e}")
        import traceback
        traceback.print_exc()
print()

print("=" * 60)
print("SECTION 3: Dependency Walker - try loading each dep individually")
print("=" * 60)
kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
LoadLibraryExW = kernel32.LoadLibraryExW
LoadLibraryExW.argtypes = [ctypes.wintypes.LPCWSTR, ctypes.wintypes.HANDLE, ctypes.wintypes.DWORD]
LoadLibraryExW.restype = ctypes.wintypes.HMODULE
LOAD_LIBRARY_AS_DATAFILE = 0x00000002

deps = [
    "python312.dll",
    "VCRUNTIME140.dll",
    "KERNEL32.dll",
    "msvcp140-5f1c5dd31916990d94181e07bc3afb32.dll",
    "api-ms-win-crt-stdio-l1-1-0.dll",
    "api-ms-win-crt-runtime-l1-1-0.dll",
    "api-ms-win-crt-heap-l1-1-0.dll",
    "api-ms-win-crt-math-l1-1-0.dll",
]

# First check which are already loaded
GetModuleHandleW = kernel32.GetModuleHandleW
GetModuleHandleW.argtypes = [ctypes.wintypes.LPCWSTR]
GetModuleHandleW.restype = ctypes.wintypes.HMODULE
GetModuleFileNameW = kernel32.GetModuleFileNameW
GetModuleFileNameW.argtypes = [ctypes.wintypes.HMODULE, ctypes.wintypes.LPWSTR, ctypes.wintypes.DWORD]
GetModuleFileNameW.restype = ctypes.wintypes.DWORD

for dep in deps:
    h = GetModuleHandleW(dep)
    if h:
        buf = ctypes.create_unicode_buffer(260)
        GetModuleFileNameW(h, buf, 260)
        print(f"  {dep:50s} LOADED at {buf.value}")
    else:
        # Try loading it
        full_path = os.path.join(libs_dir, dep)
        if os.path.exists(full_path):
            try:
                ctypes.WinDLL(full_path)
                print(f"  {dep:50s} loaded from optvl.libs OK")
            except OSError as e:
                print(f"  {dep:50s} FAILED from optvl.libs: {e}")
        else:
            try:
                ctypes.WinDLL(dep)
                print(f"  {dep:50s} loaded from system OK")
            except OSError as e:
                print(f"  {dep:50s} NOT FOUND / FAILED: {e}")
print()

print("=" * 60)
print("SECTION 4: Load .pyd as data file (bypass code validation)")
print("=" * 60)
h = LoadLibraryExW(broken_pyd, None, LOAD_LIBRARY_AS_DATAFILE)
if h:
    print(f"  LOAD_LIBRARY_AS_DATAFILE: SUCCESS (handle={h})")
    print("  -> PE structure is valid, issue is with CODE loading")
else:
    err = ctypes.get_last_error()
    print(f"  LOAD_LIBRARY_AS_DATAFILE: FAILED (error {err})")
    print("  -> PE structure itself may be corrupt")
print()

print("=" * 60)
print("SECTION 5: Check Windows version")
print("=" * 60)
ver = sys.getwindowsversion()
print(f"  Windows version: {ver.major}.{ver.minor}.{ver.build}")
print(f"  Platform: {ver.platform}")
print(f"  Service pack: {ver.service_pack}")
print()

print("Done.")