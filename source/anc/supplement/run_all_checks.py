#!/usr/bin/env python3
"""Reproduce all fifteen distinct exact certificates without changing sources.

Requires Python 3.10+ and SymPy 1.14.0 for two Part I scripts. All other
certificates use the standard library. No network requests are made.
"""
from __future__ import annotations

import difflib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


CHECKS = [
    ('part1_comparison/verify_nakayama_hom_table.py', 'part1_comparison/verify_nakayama_hom_table_output_phase4.txt'),
    ('part1_comparison/verify_nakayama_hom_table_char2.py', 'part1_comparison/verify_nakayama_hom_table_char2_output.txt'),
    ('part1_comparison/verify_a63_stable_kernel.py', 'part1_comparison/verify_a63_stable_kernel_output_phase4.txt'),
    ('part1_comparison/verify_perfect_complex_duality.py', 'part1_comparison/verify_perfect_complex_duality_output_phase4.txt'),
    ('part1_comparison/verify_ahkly_sides.py', 'part1_comparison/verify_ahkly_sides_output_phase4.txt'),
    ('part1_comparison/verify_factorisation_toys.py', 'part1_comparison/verify_factorisation_toys_output_phase4.txt'),
    *[(f'part2_fixed_kernel/scripts/{name}.py', f'part2_fixed_kernel/outputs/{name}.txt')
      for name in ('verify_dg_boundary', 'verify_nakayama_profiles', 'check_7dim',
                   'verify_formal_realisation', 'verify_group_block_obstruction')],
    ('part3_two_term/verify_nonformal.py', 'part3_two_term/verify_nonformal_output.txt'),
    ('part4_optimality/verify_connective_obstruction.py', 'part4_optimality/verification_output.txt'),
    ('part4_optimality/verify_marking_boundary.py', 'part4_optimality/marking_boundary_output.txt'),
    ('part5_revision/verify_coordinate_swap.py', 'part5_revision/verification_output.txt'),
]


def copy_materialized(source, target):
    # Read bytes explicitly so cloud-backed files are hydrated before staging.
    data = Path(source).read_bytes()
    Path(target).write_bytes(data)
    if Path(target).read_bytes() != data:
        raise OSError(f'Staging integrity check failed: {Path(source).name}')
    return target


def main() -> int:
    if sys.version_info < (3, 10):
        print('Python 3.10 or newer is required.', file=sys.stderr)
        return 2
    try:
        import sympy
    except ImportError:
        print('Install SymPy 1.14.0, or run: uv run --with sympy==1.14.0 python run_all_checks.py', file=sys.stderr)
        return 2
    print(f'Python {sys.version.split()[0]}; SymPy {sympy.__version__}', flush=True)
    root = Path(__file__).resolve().parent
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONHASHSEED='0')
    failures = 0
    with tempfile.TemporaryDirectory(prefix='bounded-defects-checks-') as tmp:
        work = Path(tmp) / 'supplement'
        shutil.copytree(root, work, copy_function=copy_materialized,
                        ignore=shutil.ignore_patterns('__pycache__', '.DS_Store', '.venv'))
        for script, expected in CHECKS:
            path = work / script
            try:
                result = subprocess.run([sys.executable, '-B', str(path)], cwd=path.parent,
                                        env=env, capture_output=True, timeout=600)
            except (OSError, subprocess.TimeoutExpired) as exc:
                print(f'{script}: FAIL ({type(exc).__name__})', flush=True)
                failures += 1
                continue
            reference = (root / expected).read_bytes()
            if result.returncode:
                print(f'{script}: FAIL (exit {result.returncode})', flush=True)
                print(result.stderr.decode(errors='replace'))
                failures += 1
            elif result.stdout != reference:
                print(f'{script}: FAIL (captured output differs)', flush=True)
                print(''.join(difflib.unified_diff(
                    reference.decode(errors='replace').splitlines(True),
                    result.stdout.decode(errors='replace').splitlines(True),
                    fromfile='captured', tofile='rerun')))
                failures += 1
            else:
                print(f'{script}: PASS (output reproduced byte-for-byte)', flush=True)
        for name in ('sphere_schur_multiplication.json', 'exterior_schur_multiplication.json'):
            rel = Path('part4_optimality') / name
            same = (work / rel).read_bytes() == (root / rel).read_bytes()
            print(f'{name}: {"PASS" if same else "FAIL"} (regenerated table)', flush=True)
            failures += not same
    if failures:
        print(f'{failures} verification(s) failed.', file=sys.stderr)
        return 1
    print('ALL 15 CERTIFICATES PASSED; ALL CAPTURED OUTPUTS AND BOTH TABLES REPRODUCED.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
