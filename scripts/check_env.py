"""
Quick environment check for non-docking + docking prerequisites.

Runs lightweight checks for:
- Python deps: RDKit import
- Binaries: obabel (Open Babel), Vina
- Network (optional): PDBe API reachability

Prints clear PASS/FAIL lines and returns non-zero on hard failures.
"""

from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path


def main() -> int:
    failures = 0

    # RDKit
    try:
        from rdkit import Chem  # type: ignore
        print("[PASS] RDKit import")
    except Exception as e:
        print(f"[FAIL] RDKit import: {e}")
        failures += 1

    # Open Babel
    obabel = shutil.which("obabel") or shutil.which("obabel.exe")
    if obabel:
        print(f"[PASS] Open Babel found: {obabel}")
    else:
        print("[WARN] Open Babel not found in PATH (required for docking)")

    # Vina binary (from env var or default name)
    vina = os.environ.get("VINA_BINARY_PATH") or shutil.which("vina") or shutil.which("vina_custom")
    if vina:
        print(f"[PASS] Vina found: {vina}")
    else:
        print("[WARN] Vina binary not found (required for docking)")

    # PDBe API (optional)
    try:
        import httpx  # type: ignore
        r = httpx.get("https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/1crn", timeout=5.0)
        if r.status_code == 200:
            print("[PASS] PDBe API reachable")
        else:
            print(f"[WARN] PDBe API status: {r.status_code}")
    except Exception as e:
        print(f"[WARN] PDBe connectivity: {e}")

    return 0 if failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

