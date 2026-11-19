#!/usr/bin/env python3
"""Run every Phase 1 certificate and check the captured output exactly.

Requires Python 3.10+ and only its standard library. No network access.
This runner verifies finite certificates; the accompanying document contains
all-degree mathematical arguments and the scope of each computation.
"""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys

NAMES = (
    "verify_dg_boundary",
    "verify_nakayama_profiles",
    "check_7dim",
    "verify_formal_realisation",
    "verify_group_block_obstruction",
)


def main() -> int:
    if sys.version_info < (3, 10):
        print("Python 3.10 or newer is required.", file=sys.stderr)
        return 2
    root = Path(__file__).resolve().parent
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PYTHONHASHSEED="0")
    failures = 0
    for name in NAMES:
        script = root / "scripts" / f"{name}.py"
        output = root / "outputs" / f"{name}.txt"
        if not script.is_file() or not output.is_file():
            print(f"{name}: FAIL (missing script or captured output)")
            failures += 1
            continue
        try:
            proc = subprocess.run(
                [sys.executable, "-B", str(script)],
                cwd=root,
                env=env,
                capture_output=True,
                check=False,
                timeout=120,
            )
        except (OSError, subprocess.TimeoutExpired) as exc:
            print(f"{name}: FAIL ({exc})")
            failures += 1
            continue
        if proc.returncode:
            print(f"{name}: FAIL (exit {proc.returncode})")
            print(proc.stderr.decode("utf-8", errors="replace"), end="")
            failures += 1
        elif proc.stdout != output.read_bytes():
            print(f"{name}: FAIL (captured output differs)")
            print(proc.stdout.decode("utf-8", errors="replace"), end="")
            failures += 1
        else:
            print(f"{name}: PASS (output reproduced byte-for-byte)")
    if failures:
        print(f"{failures} certificate(s) failed.", file=sys.stderr)
        return 1
    print("ALL FIVE CERTIFICATES PASSED AND REPRODUCED THEIR CAPTURED OUTPUTS.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
