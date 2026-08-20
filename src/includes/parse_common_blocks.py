#!/usr/bin/env python3
"""
Parse a Fortran common-block include file and report the memory footprint
of every variable and every COMMON block.

Parameter constants (NVMAX, NSMAX, …) are read automatically from
INCLUDE'd files (ADIMEN.INC, AINDEX.INC, etc.) and from PARAMETER
statements in the main file itself.

Supports a multiplier per file — use FILE:COUNT syntax to indicate how
many times a file's storage is replicated.

Usage:
    python parse_common_blocks.py AVL.INC
    python parse_common_blocks.py AVL.INC:30 AVL_surf.INC:4
    python parse_common_blocks.py AVL_oml.INC:4 AVL_surf.INC:4 AVL.INC:30 AVL_ad_seeds.INC:14
    python parse_common_blocks.py AVL.INC -I /path/to/includes
    python parse_common_blocks.py AVL.INC -p ADIMEN.INC AINDEX.INC
    python parse_common_blocks.py AVL.INC --avl-real 4
    python parse_common_blocks.py AVL.INC --set NVMAX 10000 --set NOBMAX 5000
    
    Here is what I use currently to get the total common block allocation with llvm-flang
    There seems to be some issue with llvm-flang that causes all the commonblocks in each file adds data independently
    python parse_common_blocks.py AVL_oml.INC:4 AVL_surf.INC:4 AVL.INC:30 AVL_ad_seeds.INC:14 AVL_surf_ad_seeds.INC:2
"""

import argparse
import re
import sys
import math
from collections import OrderedDict
from pathlib import Path

# ============================================================
# Fallback defaults — used only when an include file cannot be
# found or a parameter is referenced but never defined.
# ============================================================
FALLBACK_PARAMS: dict[str, int] = {
    "NVMAX":   5000,
    "NSMAX":   500,
    "NSECMAX": 301,
    "NFMAX":   100,
    "NLMAX":   500,
    "NBMAX":   20,

    "NUMAX":   6,
    "NDMAX":   30,
    "NGMAX":   21,

    "NRMAX":   25,
    "NTMAX":   503,

    "NOBMAX":  1,

    "ICONX":   20,
    "IBX":     200,

    "IVTOT":   20,
    "ICTOT":   20,
    "IPTOT":   20,
    "JETOT":   12,
    "KPTOT":   20,
    "KPVTOT":  10,
}


# ============================================================
# Type → bytes mapping
# ============================================================

def make_type_bytes(avl_real_bytes: int) -> dict[str, int]:
    return {
        "integer":          4,
        "real":             avl_real_bytes,
        "real*4":           4,
        "real*8":           8,
        "real(kind=avl_real)": avl_real_bytes,
        "double precision": 8,
        "logical":          4,
        "complex":          2 * avl_real_bytes,
        "complex*8":        8,
        "complex*16":       16,
        "complex(kind=avl_real)": 2 * avl_real_bytes,
    }


# ============================================================
# Helpers
# ============================================================

def eval_dim(expr: str, params: dict[str, int]) -> int:
    """Safely evaluate a dimension expression using *params*."""
    expr = expr.strip()
    safe = re.sub(r"[A-Za-z_]\w*", lambda m: str(params.get(m.group(), m.group())), expr)
    try:
        return int(eval(safe, {"__builtins__": {}}))
    except Exception as e:
        print(f"  *** WARNING: cannot evaluate dimension '{expr}' -> '{safe}': {e}",
              file=sys.stderr)
        return 0


def human_bytes(n: int) -> str:
    """Return a human-readable byte string."""
    n = int(n)
    if n == 0:
        return "0 B"
    units = ["B", "KiB", "MiB", "GiB", "TiB"]
    idx = min(int(math.log2(n) // 10), len(units) - 1) if n > 0 else 0
    value = n / (1 << (10 * idx))
    if idx == 0:
        return f"{n} B"
    return f"{value:,.2f} {units[idx]}"


# ============================================================
# Fortran source reader — joins continuation lines
# ============================================================

def read_fortran_lines(path: str) -> list[str]:
    """Read a fixed-form Fortran file, join continuation lines,
    strip comments, and return a list of logical lines."""
    raw = Path(path).read_text(errors="replace").splitlines()
    logical: list[str] = []

    for line in raw:
        # Completely blank or comment line (C/c/*/! in column 1)
        if not line or line[0] in ("C", "c", "!", "*"):
            continue

        # Inline comment after '!'
        bang = line.find("!")
        if bang > 0:
            line = line[:bang]

        # Fixed-form continuation: non-blank, non-zero char in column 6
        if len(line) > 5 and line[0] == " " and line[5] not in (" ", "0", "\t"):
            if logical:
                logical[-1] += " " + line[6:].strip()
            continue

        # Tab-form continuation
        if line and line[0] == "\t" and len(line) > 1 and line[1].isdigit():
            if logical:
                logical[-1] += " " + line[2:].strip()
            continue

        logical.append(line.rstrip())

    return logical


# ============================================================
# Parse PARAMETER statements from Fortran source
# ============================================================

_RE_PARAMETER = re.compile(
    r"^\s*PARAMETER\s*\((.+)\)\s*$", re.IGNORECASE
)

def parse_parameters_from_lines(lines: list[str], params: dict[str, int]):
    for line in lines:
        m = _RE_PARAMETER.match(line)
        if not m:
            continue

        body = m.group(1)
        depth = 0
        parts: list[str] = []
        cur: list[str] = []
        for ch in body:
            if ch == "(":
                depth += 1; cur.append(ch)
            elif ch == ")":
                depth -= 1; cur.append(ch)
            elif ch == "," and depth == 0:
                parts.append("".join(cur).strip()); cur = []
            else:
                cur.append(ch)
        parts.append("".join(cur).strip())

        for part in parts:
            if "=" not in part:
                continue
            name, expr = part.split("=", 1)
            name = name.strip().upper()
            expr = expr.strip()
            val = eval_dim(expr, params)
            if val != 0 or expr.strip() == "0":
                params[name] = val


def parse_parameters_from_file(path: str, params: dict[str, int]):
    lines = read_fortran_lines(path)
    parse_parameters_from_lines(lines, params)


# ============================================================
# Resolve INCLUDE directives
# ============================================================

_RE_INCLUDE = re.compile(
    r"""^\s*INCLUDE\s+['"]([^'"]+)['"]\s*$""", re.IGNORECASE
)

def find_include_file(name: str, search_dirs: list[Path]) -> Path:
    for d in search_dirs:
        candidate = d / name
        if candidate.is_file():
            return candidate
    return None


def resolve_includes(
    main_path: str,
    extra_search_dirs: list[str] = None,
) -> list[str]:
    main_dir = Path(main_path).resolve().parent
    search_dirs = [main_dir]
    if extra_search_dirs:
        search_dirs.extend(Path(d).resolve() for d in extra_search_dirs)

    raw = Path(main_path).read_text(errors="replace").splitlines()
    found: list[str] = []
    not_found: list[str] = []

    for line in raw:
        stripped = line.lstrip()
        if stripped and stripped[0] in ("C", "c", "!", "*"):
            continue
        m = _RE_INCLUDE.match(line)
        if m:
            inc_name = m.group(1)
            inc_path = find_include_file(inc_name, search_dirs)
            if inc_path:
                found.append(str(inc_path))
            else:
                not_found.append(inc_name)

    return found, not_found


# ============================================================
# Parse type declarations
# ============================================================

_RE_TYPE_DECL = re.compile(
    r"""
    ^\s*
    (?:
        (INTEGER|LOGICAL)
      | (REAL|COMPLEX)\s*\(\s*kind\s*=\s*avl_real\s*\)
      | (REAL\*\d+|COMPLEX\*\d+|DOUBLE\s+PRECISION)
      | (REAL|COMPLEX|INTEGER|LOGICAL)
      | (CHARACTER)\*(\d+)
    )
    \s+
    (.+)
    """,
    re.IGNORECASE | re.VERBOSE,
)

def parse_type_declarations(lines: list[str], type_bytes: dict[str, int],
                            avl_real_bytes: int) -> dict[str, tuple[str, int]]:
    var_types: dict[str, tuple[str, int]] = {}

    for line in lines:
        m = _RE_TYPE_DECL.match(line)
        if not m:
            continue

        if m.group(1):
            tkey = m.group(1).upper().strip()
        elif m.group(2):
            tkey = m.group(2).lower() + "(kind=avl_real)"
        elif m.group(3):
            tkey = re.sub(r"\s+", " ", m.group(3)).lower().strip()
        elif m.group(4):
            tkey = m.group(4).lower().strip()
        elif m.group(5):
            tkey = f"character*{m.group(6)}"
        else:
            continue

        bpe = type_bytes.get(tkey.lower())
        if bpe is None:
            cm = re.match(r"character\*(\d+)", tkey, re.IGNORECASE)
            if cm:
                bpe = int(cm.group(1))
            else:
                bpe = avl_real_bytes

        names_str = m.group(7) if m.group(7) else ""
        if not names_str:
            continue
        depth = 0
        parts: list[str] = []
        cur: list[str] = []
        for ch in names_str:
            if ch == "(":
                depth += 1; cur.append(ch)
            elif ch == ")":
                depth -= 1; cur.append(ch)
            elif ch == "," and depth == 0:
                parts.append("".join(cur).strip()); cur = []
            else:
                cur.append(ch)
        parts.append("".join(cur).strip())

        for part in parts:
            name = re.split(r"[(\s]", part)[0].strip().upper()
            if name:
                var_types[name] = (tkey, bpe)

    return var_types


# ============================================================
# Parse CHARACTER declarations (old-style)
# ============================================================

_RE_CHAR_DECL = re.compile(
    r"^\s*CHARACTER\*(\d+)\s+(.+)", re.IGNORECASE
)

def parse_character_decls(lines: list[str], var_types: dict[str, tuple[str, int]]):
    for line in lines:
        m = _RE_CHAR_DECL.match(line)
        if not m:
            continue
        charlen = int(m.group(1))
        names_str = m.group(2)
        depth = 0
        parts: list[str] = []
        cur: list[str] = []
        for ch in names_str:
            if ch == "(":
                depth += 1; cur.append(ch)
            elif ch == ")":
                depth -= 1; cur.append(ch)
            elif ch == "," and depth == 0:
                parts.append("".join(cur).strip()); cur = []
            else:
                cur.append(ch)
        parts.append("".join(cur).strip())
        for part in parts:
            name = re.split(r"[(\s]", part)[0].strip().upper()
            if name:
                var_types[name] = (f"character*{charlen}", charlen)


# ============================================================
# Parse COMMON blocks
# ============================================================

_RE_COMMON = re.compile(r"^\s*COMMON\s*/\s*(\w+)\s*/\s*(.+)", re.IGNORECASE)

def parse_common_blocks(lines: list[str]) -> OrderedDict[str, list[tuple[str, str]]]:
    blocks: OrderedDict[str, list[tuple[str, str]]] = OrderedDict()
    current_block: str = None
    current_tail = ""

    for line in lines:
        m = _RE_COMMON.match(line)
        if m:
            if current_block is not None:
                _parse_var_list(current_tail, blocks[current_block])
            current_block = m.group(1).upper()
            if current_block not in blocks:
                blocks[current_block] = []
            current_tail = m.group(2)
            continue

    if current_block is not None:
        _parse_var_list(current_tail, blocks[current_block])

    return blocks


def _parse_var_list(text: str, out: list[tuple[str, str]]):
    text = re.sub(r"!.*", "", text)
    depth = 0
    parts: list[str] = []
    cur: list[str] = []
    for ch in text:
        if ch == "(":
            depth += 1; cur.append(ch)
        elif ch == ")":
            depth -= 1; cur.append(ch)
        elif ch == "," and depth == 0:
            parts.append("".join(cur).strip()); cur = []
        else:
            cur.append(ch)
    tail = "".join(cur).strip()
    if tail:
        parts.append(tail)

    for part in parts:
        part = part.strip()
        if not part:
            continue
        pm = re.match(r"(\w+)\s*\(([^)]*)\)", part)
        if pm:
            out.append((pm.group(1).upper(), pm.group(2)))
        else:
            name = part.strip().upper()
            if name and re.match(r"[A-Z_]\w*$", name):
                out.append((name, ""))


# ============================================================
# Compute sizes
# ============================================================

def compute_sizes(
    blocks: OrderedDict[str, list[tuple[str, str]]],
    var_types: dict[str, tuple[str, int]],
    params: dict[str, int],
    avl_real_bytes: int,
) -> OrderedDict[str, list[dict]]:
    result: OrderedDict[str, list[dict]] = OrderedDict()

    for bname, varlist in blocks.items():
        entries: list[dict] = []
        for vname, dims_str in varlist:
            tinfo = var_types.get(vname)
            if tinfo:
                tkey, bpe = tinfo
            else:
                if bname.endswith("_I"):
                    tkey, bpe = "integer", 4
                elif bname.endswith("_L"):
                    tkey, bpe = "logical", 4
                else:
                    tkey, bpe = "real(kind=avl_real)", avl_real_bytes

            if dims_str:
                dim_parts = []
                depth = 0
                cur: list[str] = []
                for ch in dims_str:
                    if ch == "(":
                        depth += 1; cur.append(ch)
                    elif ch == ")":
                        depth -= 1; cur.append(ch)
                    elif ch == "," and depth == 0:
                        dim_parts.append("".join(cur).strip()); cur = []
                    else:
                        cur.append(ch)
                dim_parts.append("".join(cur).strip())

                nelems = 1
                dim_vals: list[int] = []
                for d in dim_parts:
                    v = eval_dim(d, params)
                    dim_vals.append(v)
                    nelems *= v
            else:
                nelems = 1
                dim_vals = []

            total_bytes = nelems * bpe
            entries.append({
                "name": vname,
                "type": tkey,
                "bpe": bpe,
                "dims_str": dims_str,
                "dim_vals": dim_vals,
                "nelems": nelems,
                "bytes": total_bytes,
            })

        result[bname] = entries

    return result


# ============================================================
# Pretty-print — weighted report
# ============================================================

def print_weighted_report(
    file_results: list[dict],
    params: dict[str, int],
):
    """Print a combined report across all files, weighted by their counts.

    file_results is a list of dicts:
        {
            "path": str,
            "count": int,
            "result": OrderedDict[block_name -> [entries]],
        }
    """

    # ── Per-file summary ──
    print(f"\n{'='*90}")
    print(f"  PER-FILE SUMMARY")
    print(f"{'='*90}")
    print(f"  {'File':<30s} {'Count':>6s} {'Single copy':>14s} {'Total (× count)':>18s}")
    print(f"  {'-'*30} {'-'*6} {'-'*14} {'-'*18}")

    overall_total = 0
    for fr in file_results:
        single = sum(e["bytes"] for entries in fr["result"].values() for e in entries)
        weighted = single * fr["count"]
        overall_total += weighted
        print(f"  {Path(fr['path']).name:<30s} {fr['count']:>6d} {human_bytes(single):>14s} {human_bytes(weighted):>18s}")

    print(f"  {'-'*30} {'-'*6} {'-'*14} {'-'*18}")
    print(f"  {'GRAND TOTAL':<30s} {'':>6s} {'':>14s} {human_bytes(overall_total):>18s}")

    # ── Per-file block breakdown ──
    for fr in file_results:
        fname = Path(fr["path"]).name
        count = fr["count"]
        result = fr["result"]
        file_single = sum(e["bytes"] for entries in result.values() for e in entries)
        file_total = file_single * count

        print(f"\n{'='*90}")
        print(f"  COMMON BLOCKS IN: {fname}  (×{count})    single={human_bytes(file_single)}  total={human_bytes(file_total)}")
        print(f"{'='*90}")
        print(f"  {'Block':<30s} {'Single copy':>14s} {'× {:<4d}'.format(count):>12s} {'% of grand':>10s}")
        print(f"  {'-'*30} {'-'*14} {'-'*12} {'-'*10}")

        block_list = []
        for bname, entries in result.items():
            bsize = sum(e["bytes"] for e in entries)
            block_list.append((bname, bsize))

        block_list.sort(key=lambda x: -x[1])
        for bname, bsize in block_list:
            weighted = bsize * count
            pct = 100.0 * weighted / overall_total if overall_total else 0
            print(f"  {bname:<30s} {human_bytes(bsize):>14s} {human_bytes(weighted):>12s} {pct:>9.1f}%")

        print(f"  {'-'*30} {'-'*14} {'-'*12} {'-'*10}")
        pct_file = 100.0 * file_total / overall_total if overall_total else 0
        print(f"  {'FILE TOTAL':<30s} {human_bytes(file_single):>14s} {human_bytes(file_total):>12s} {pct_file:>9.1f}%")

    # ── Weighted COMMON block ranking (across all files) ──
    # A block name may appear in multiple files; aggregate weighted sizes.
    block_weighted: dict[str, int] = {}
    block_files: dict[str, list[str]] = {}
    for fr in file_results:
        fname = Path(fr["path"]).name
        count = fr["count"]
        for bname, entries in fr["result"].items():
            bsize = sum(e["bytes"] for e in entries)
            weighted = bsize * count
            block_weighted[bname] = block_weighted.get(bname, 0) + weighted
            block_files.setdefault(bname, []).append(f"{fname}×{count}")

    sorted_blocks = sorted(block_weighted.items(), key=lambda x: -x[1])

    print(f"\n{'='*90}")
    print(f"  WEIGHTED COMMON BLOCK RANKING (all files combined)")
    print(f"{'='*90}")
    print(f"  {'#':<4s} {'Block':<24s} {'Weighted size':>14s} {'% of total':>10s} {'Sources':<30s}")
    print(f"  {'-'*4} {'-'*24} {'-'*14} {'-'*10} {'-'*30}")

    for rank, (bname, wsize) in enumerate(sorted_blocks[:40], 1):
        pct = 100.0 * wsize / overall_total if overall_total else 0
        sources = ", ".join(block_files[bname])
        print(f"  {rank:<4d} {bname:<24s} {human_bytes(wsize):>14s} {pct:>9.1f}% {sources:<30s}")

    # ── Top variables by weighted size ──
    # Collect (file, count, block, var_entry) tuples
    all_vars: list[tuple[str, int, str, dict]] = []
    for fr in file_results:
        fname = Path(fr["path"]).name
        count = fr["count"]
        for bname, entries in fr["result"].items():
            for e in entries:
                all_vars.append((fname, count, bname, e))

    # Sort by weighted bytes (size × count)
    all_vars.sort(key=lambda x: -(x[3]["bytes"] * x[1]))
    top_n = all_vars[:50]

    if top_n:
        print(f"\n{'='*110}")
        print(f"  TOP {len(top_n)} LARGEST VARIABLES (by weighted size = single size × count)")
        print(f"{'='*110}")
        print(f"  {'#':<4s} {'Variable':<20s} {'Block':<18s} {'File':<20s} {'Cnt':>4s} {'Type':<22s} "
              f"{'Single size':>12s} {'Weighted':>14s} {'% total':>8s}")
        print(f"  {'-'*4} {'-'*20} {'-'*18} {'-'*20} {'-'*4} {'-'*22} "
              f"{'-'*12} {'-'*14} {'-'*8}")

        for rank, (fname, count, bname, e) in enumerate(top_n, 1):
            weighted = e["bytes"] * count
            pct = 100.0 * weighted / overall_total if overall_total else 0
            print(f"  {rank:<4d} {e['name']:<20s} {bname:<18s} {fname:<20s} {count:>4d} {e['type']:<22s} "
                  f"{human_bytes(e['bytes']):>12s} {human_bytes(weighted):>14s} {pct:>7.1f}%")

        print(f"  {'-'*4} {'-'*20} {'-'*18} {'-'*20} {'-'*4} {'-'*22} "
              f"{'-'*12} {'-'*14} {'-'*8}")
        top_total = sum(x[3]["bytes"] * x[1] for x in top_n)
        top_pct = 100.0 * top_total / overall_total if overall_total else 0
        print(f"  {'':4s} {'Top ' + str(len(top_n)) + ' total':<20s} {'':18s} {'':20s} {'':>4s} {'':22s} "
              f"{'':>12s} {human_bytes(top_total):>14s} {top_pct:>7.1f}%")
        
    print(f"  {'-'*80}")
    print(f"  {'GRAND TOTAL':<30s} {'':>6s} {'':>14s} {human_bytes(overall_total):>18s}")
    print(f"  {'-'*80}")


    print()


# ============================================================
# Main
# ============================================================

def parse_file_spec(spec: str) -> tuple[str, int]:
    """Parse 'filename:count' or just 'filename' (count defaults to 1)."""
    # Handle Windows paths with drive letters like C:\...:4
    # Try splitting from the right
    parts = spec.rsplit(":", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0], int(parts[1])
    return spec, 1


def main():
    ap = argparse.ArgumentParser(
        description="Report weighted memory footprint of Fortran COMMON blocks.  "
                    "Use FILE:COUNT syntax to specify how many times each file's "
                    "storage is replicated (e.g. AVL.INC:30).  "
                    "The report shows totals weighted by count.",
        epilog="Example:\n"
               "  python parse_common_blocks.py AVL_oml.INC:4 AVL_surf.INC:4 "
               "AVL.INC:30 AVL_ad_seeds.INC:14\n",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "include_files",
        nargs="+",
        help="Fortran include file(s) with optional :COUNT suffix "
             "(e.g. AVL.INC:30).  Count defaults to 1.",
    )
    ap.add_argument(
        "-p", "--params",
        nargs="+", metavar="FILE",
        help="Fortran file(s) to read PARAMETER values from.",
    )
    ap.add_argument(
        "-I", "--include-dir",
        action="append", default=[], metavar="DIR",
        help="Additional directory to search for INCLUDE'd files.",
    )
    ap.add_argument(
        "--avl-real",
        type=int, default=8, choices=[4, 8],
        help="Byte width of avl_real kind (default: 8 = double precision).",
    )
    ap.add_argument(
        "-d", "--details",
        action="store_true",
        help="Print per-variable breakdown for each COMMON block per file.",
    )
    ap.add_argument(
        "--set",
        nargs=2, action="append", default=[], metavar=("NAME", "VALUE"),
        help="Override or add a parameter value, e.g. --set NVMAX 10000.",
    )
    args = ap.parse_args()

    avl_real_bytes: int = args.avl_real
    type_bytes = make_type_bytes(avl_real_bytes)

    # Parse file:count specs
    file_specs: list[tuple[str, int]] = []
    for spec in args.include_files:
        path, count = parse_file_spec(spec)
        file_specs.append((path, count))

    input_paths = [p for p, _ in file_specs]

    for p in input_paths:
        if not Path(p).is_file():
            print(f"Error: file not found: {p}", file=sys.stderr)
            sys.exit(1)

    # ----------------------------------------------------------
    # 1. Build a shared parameter table
    # ----------------------------------------------------------
    params: dict[str, int] = dict(FALLBACK_PARAMS)

    param_files_found: list[str] = []
    param_files_missing: list[str] = []

    if args.params:
        search_dirs: list[Path] = []
        for ip in input_paths:
            d = Path(ip).resolve().parent
            if d not in search_dirs:
                search_dirs.append(d)
        search_dirs.extend(Path(d).resolve() for d in args.include_dir)

        for pf in args.params:
            p = Path(pf)
            if p.is_file():
                param_files_found.append(str(p))
            else:
                found_path = find_include_file(p.name, search_dirs)
                if found_path:
                    param_files_found.append(str(found_path))
                else:
                    param_files_missing.append(pf)
    else:
        seen: set[str] = set()
        for ip in input_paths:
            found, not_found = resolve_includes(ip, args.include_dir)
            for f in found:
                if f not in seen:
                    param_files_found.append(f)
                    seen.add(f)
            for nf in not_found:
                if nf not in seen:
                    param_files_missing.append(nf)
                    seen.add(nf)

    for pf in param_files_found:
        parse_parameters_from_file(pf, params)

    all_lines: dict[str, list[str]] = {}
    for ip in input_paths:
        lines = read_fortran_lines(ip)
        all_lines[ip] = lines
        parse_parameters_from_lines(lines, params)

    for name, value in args.set:
        params[name.upper()] = int(value)

    # ----------------------------------------------------------
    # 2. Process each file independently, keeping its count
    # ----------------------------------------------------------
    file_results: list[dict] = []

    for ip, count in file_specs:
        lines = all_lines[ip]

        # Types
        var_types: dict[str, tuple[str, int]] = {}
        parse_character_decls(lines, var_types)
        var_types.update(parse_type_declarations(lines, type_bytes, avl_real_bytes))

        # COMMON blocks
        blocks = parse_common_blocks(lines)
        result = compute_sizes(blocks, var_types, params, avl_real_bytes)

        file_results.append({
            "path": ip,
            "count": count,
            "result": result,
        })

    # ----------------------------------------------------------
    # 3. Print header
    # ----------------------------------------------------------
    print(f"\nFortran Common Block Weighted Memory Report")
    print(f"{'='*60}")
    print(f"  {'File':<30s} {'Count':>6s}")
    print(f"  {'-'*30} {'-'*6}")
    for ip, count in file_specs:
        print(f"  {Path(ip).name:<30s} {count:>6d}")
    print(f"avl_real = {avl_real_bytes} bytes")

    if param_files_found:
        print(f"\nParameter files read:")
        for pf in param_files_found:
            print(f"  {pf}")

    if param_files_missing:
        print(f"\nWARNING — include files not found (using fallback values):")
        for pf in param_files_missing:
            print(f"  {pf}")

    print(f"\nParameter values used:")
    for k, v in sorted(params.items()):
        print(f"  {k:>12s} = {v:>8,d}")

    # ----------------------------------------------------------
    # 4. Per-file detailed breakdown (if requested)
    # ----------------------------------------------------------
    if args.details:
        for fr in file_results:
            fname = Path(fr["path"]).name
            count = fr["count"]
            result = fr["result"]

            for bname, entries in result.items():
                block_total = sum(e["bytes"] for e in entries)
                print(f"\n{'='*100}")
                print(f"  COMMON /{bname}/  from {fname} ×{count}    "
                      f"(single={human_bytes(block_total)}, weighted={human_bytes(block_total * count)})")
                print(f"{'='*100}")
                print(f"  {'Variable':<20s} {'Type':<26s} {'Dimensions':<30s} "
                      f"{'Elements':>12s} {'Single':>14s} {'Weighted':>14s}")
                print(f"  {'-'*20} {'-'*26} {'-'*30} {'-'*12} {'-'*14} {'-'*14}")

                for e in entries:
                    dims_display = ""
                    if e["dim_vals"]:
                        dims_display = "(" + ", ".join(str(v) for v in e["dim_vals"]) + ")"
                        if e["dims_str"] and len(dims_display) < 20:
                            dims_display += f"  [{e['dims_str']}]"

                    weighted = e["bytes"] * count
                    print(f"  {e['name']:<20s} {e['type']:<26s} {dims_display:<30s} "
                          f"{e['nelems']:>12,d} {human_bytes(e['bytes']):>14s} {human_bytes(weighted):>14s}")

                print(f"  {'':->116s}")
                print(f"  {'Block total':>90s}: {human_bytes(block_total):>14s} → {human_bytes(block_total * count)}")

    # ----------------------------------------------------------
    # 5. Weighted summary report
    # ----------------------------------------------------------
    print_weighted_report(file_results, params)


if __name__ == "__main__":
    main()