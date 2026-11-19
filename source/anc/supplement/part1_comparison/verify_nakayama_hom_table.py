#!/usr/bin/env python3
"""Verify the cone-Hom table for A_d = k C_{2d}/rad^d, d=3,4,5.

Right-module convention:
  P_i=e_i A, C_i=[P_{i+1} --alpha_i--> P_i] in degrees -1,0.
  Hom_A(P_a,P_b)=e_b A e_a, with basis the unique path b -> a of
  length (a-b mod 2d) when that length is < d.

Mapping complex:
  C^{-1}=Hom(P_i,P_{j+1})
  C^0 =Hom(P_{i+1},P_{j+1}) direct-sum Hom(P_i,P_j)
  C^1 =Hom(P_{i+1},P_j)
  delta^{-1}(f)=(f alpha_i, alpha_j f)
  delta^0(a,b)=alpha_j a - b alpha_i.

The script constructs the matrices over Q, checks delta^0 delta^{-1}=0,
and compares the resulting cohomology with the claimed exceptional offsets
(0,0), (1,1), (d,-1), (d+1,0).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple
import sympy as sp


def path_exists(start: int, end: int, n: int, d: int) -> bool:
    """Unique oriented path start -> end survives iff its length is < d."""
    return ((end - start) % n) < d


def path_length(start: int, end: int, n: int) -> int:
    return (end - start) % n


@dataclass(frozen=True)
class Row:
    d: int
    m: int
    dims: Tuple[int, int, int]
    dminus: sp.Matrix
    dzero: sp.Matrix
    ranks: Tuple[int, int]
    h: Tuple[int, int, int]


def mapping_row(d: int, m: int) -> Row:
    n = 2 * d
    i = 0
    j = m % n

    # Basis dimensions of the three mapping-complex terms.
    cminus = int(path_exists(j + 1, i, n, d))
    c0a = int(path_exists(j + 1, i + 1, n, d))
    c0b = int(path_exists(j, i, n, d))
    c0 = c0a + c0b
    cplus = int(path_exists(j, i + 1, n, d))

    # delta^{-1}: each nonzero composition is the corresponding path with + sign.
    dm = sp.zeros(c0, cminus)
    if cminus:
        f_len = path_length(j + 1, i, n)
        row = 0
        if c0a:
            # f alpha_i: j+1 -> i -> i+1 (path length f_len+1)
            if f_len + 1 < d:
                dm[row, 0] = 1
            row += 1
        if c0b:
            # alpha_j f: j -> j+1 -> i (path length f_len+1)
            if f_len + 1 < d:
                dm[row, 0] = 1

    # delta^0(a,b)=alpha_j a - b alpha_i.
    dz = sp.zeros(cplus, c0)
    if cplus:
        target_len = path_length(j, i + 1, n)
        col = 0
        if c0a:
            a_len = path_length(j + 1, i + 1, n)
            if a_len + 1 < d and a_len + 1 == target_len:
                dz[0, col] = 1
            col += 1
        if c0b:
            b_len = path_length(j, i, n)
            if b_len + 1 < d and b_len + 1 == target_len:
                dz[0, col] = -1

    if dz * dm != sp.zeros(cplus, cminus):
        raise AssertionError(f'd^2 != 0 for d={d}, m={m}: {dz*dm}')

    rminus = dm.rank()
    rzero = dz.rank()
    hminus = cminus - rminus
    hzero = c0 - rminus - rzero
    hplus = cplus - rzero
    return Row(d, m, (cminus, c0, cplus), dm, dz,
               (rminus, rzero), (hminus, hzero, hplus))


def expected_h(d: int, m: int) -> Tuple[int, int, int]:
    if m == 0:
        return (0, 1, 0)
    if m == 1:
        return (0, 0, 1)
    if m == d:
        return (1, 0, 0)
    if m == d + 1:
        return (0, 1, 0)
    return (0, 0, 0)


def fmt_matrix(M: sp.Matrix) -> str:
    if M.rows == 0 or M.cols == 0:
        return f'{M.rows}x{M.cols}:[]'
    return f'{M.rows}x{M.cols}:' + str(M.tolist()).replace(' ', '')


def main() -> None:
    all_ok = True
    for d in (3, 4, 5):
        print(f'=== d={d}, n={2*d} ===')
        print('m | dims(C-1,C0,C1) | d-1 | d0 | ranks | H(-1,0,1) | expected | status')
        for m in range(2 * d):
            row = mapping_row(d, m)
            exp = expected_h(d, m)
            ok = row.h == exp
            all_ok &= ok
            print(
                f'{m:2d} | {row.dims!s:>14} | {fmt_matrix(row.dminus):>14} | '
                f'{fmt_matrix(row.dzero):>14} | {row.ranks!s:>6} | '
                f'{row.h!s:>9} | {exp!s:>9} | {"OK" if ok else "FAIL"}'
            )
        print()
    if not all_ok:
        raise SystemExit('Verification failed.')
    print('ALL CHECKS PASSED for d=3,4,5.')


if __name__ == '__main__':
    main()
