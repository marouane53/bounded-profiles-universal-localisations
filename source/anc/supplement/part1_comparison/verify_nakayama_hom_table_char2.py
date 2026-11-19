#!/usr/bin/env python3
"""Verify the Nakayama cone-Hom table over GF(2) for d=3,4,5.

This dependency-free companion to ``verify_nakayama_hom_table.py`` constructs
the same three-term mapping complexes, reduces every coefficient modulo 2,
checks that the two differentials compose to zero, and computes cohomology
dimensions using exact Gaussian elimination over GF(2).

In characteristic two the row ``(1,-1)`` becomes ``(1,1)``.  Its rank is still
one, and its product with the column ``(1,1)^t`` is zero because 1+1=0.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Tuple


@dataclass(frozen=True)
class MatrixGF2:
    """A small matrix whose entries are normalized to 0 or 1."""

    rows: int
    cols: int
    entries: Tuple[Tuple[int, ...], ...]

    @classmethod
    def zeros(cls, rows: int, cols: int) -> "MatrixGF2":
        return cls(rows, cols, tuple(tuple(0 for _ in range(cols)) for _ in range(rows)))

    @classmethod
    def from_rows(
        cls, rows: int, cols: int, entries: Iterable[Iterable[int]]
    ) -> "MatrixGF2":
        normalized = tuple(tuple(value % 2 for value in row) for row in entries)
        if len(normalized) != rows or any(len(row) != cols for row in normalized):
            raise ValueError("matrix data do not match the declared shape")
        return cls(rows, cols, normalized)

    def multiply(self, other: "MatrixGF2") -> "MatrixGF2":
        if self.cols != other.rows:
            raise ValueError("incompatible matrix shapes")
        values = [
            [
                sum(self.entries[i][k] * other.entries[k][j] for k in range(self.cols))
                % 2
                for j in range(other.cols)
            ]
            for i in range(self.rows)
        ]
        return MatrixGF2.from_rows(self.rows, other.cols, values)

    def rank(self) -> int:
        work = [list(row) for row in self.entries]
        pivot_row = 0
        for col in range(self.cols):
            pivot = next(
                (row for row in range(pivot_row, self.rows) if work[row][col]),
                None,
            )
            if pivot is None:
                continue
            work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
            for row in range(self.rows):
                if row != pivot_row and work[row][col]:
                    work[row] = [a ^ b for a, b in zip(work[row], work[pivot_row])]
            pivot_row += 1
            if pivot_row == self.rows:
                break
        return pivot_row

    def is_zero(self) -> bool:
        return all(value == 0 for row in self.entries for value in row)

    def display(self) -> str:
        if self.rows == 0 or self.cols == 0:
            return f"{self.rows}x{self.cols}:[]"
        compact = str([list(row) for row in self.entries]).replace(" ", "")
        return f"{self.rows}x{self.cols}:{compact}"


def path_exists(start: int, end: int, n: int, d: int) -> bool:
    """The unique oriented path start -> end survives exactly below length d."""

    return ((end - start) % n) < d


def path_length(start: int, end: int, n: int) -> int:
    return (end - start) % n


def expected_h(d: int, offset: int) -> Tuple[int, int, int]:
    exceptional = {
        0: (0, 1, 0),
        1: (0, 0, 1),
        d: (1, 0, 0),
        d + 1: (0, 1, 0),
    }
    return exceptional.get(offset, (0, 0, 0))


def mapping_row(
    d: int, offset: int
) -> Tuple[
    Tuple[int, int, int],
    MatrixGF2,
    MatrixGF2,
    Tuple[int, int],
    Tuple[int, int, int],
]:
    n = 2 * d
    i = 0
    j = offset % n

    cminus = int(path_exists(j + 1, i, n, d))
    c0a = int(path_exists(j + 1, i + 1, n, d))
    c0b = int(path_exists(j, i, n, d))
    c0 = c0a + c0b
    cplus = int(path_exists(j, i + 1, n, d))

    dminus_entries = [[0 for _ in range(cminus)] for _ in range(c0)]
    if cminus:
        f_length = path_length(j + 1, i, n)
        row = 0
        if c0a:
            if f_length + 1 < d:
                dminus_entries[row][0] = 1
            row += 1
        if c0b and f_length + 1 < d:
            dminus_entries[row][0] = 1
    dminus = MatrixGF2.from_rows(c0, cminus, dminus_entries)

    dzero_entries = [[0 for _ in range(c0)] for _ in range(cplus)]
    if cplus:
        target_length = path_length(j, i + 1, n)
        col = 0
        if c0a:
            a_length = path_length(j + 1, i + 1, n)
            if a_length + 1 < d and a_length + 1 == target_length:
                dzero_entries[0][col] = 1
            col += 1
        if c0b:
            b_length = path_length(j, i, n)
            if b_length + 1 < d and b_length + 1 == target_length:
                # The rational coefficient is -1; modulo 2 this is +1.
                dzero_entries[0][col] = (-1) % 2
    dzero = MatrixGF2.from_rows(cplus, c0, dzero_entries)

    composite = dzero.multiply(dminus)
    if not composite.is_zero():
        raise AssertionError(f"d^2 != 0 for d={d}, offset={offset}: {composite}")

    rank_minus = dminus.rank()
    rank_zero = dzero.rank()
    cohomology = (
        cminus - rank_minus,
        c0 - rank_minus - rank_zero,
        cplus - rank_zero,
    )
    return (
        (cminus, c0, cplus),
        dminus,
        dzero,
        (rank_minus, rank_zero),
        cohomology,
    )


def main() -> None:
    print("Nakayama cone-Hom characteristic-two certificate (GF(2))")
    print("coefficient reduction -1 = +1: PASS")
    all_ok = True
    for d in (3, 4, 5):
        print(f"=== d={d}, n={2 * d} ===")
        print("m | dims(C-1,C0,C1) | d-1 | d0 | ranks | H(-1,0,1) | expected | status")
        for offset in range(2 * d):
            dims, dminus, dzero, ranks, cohomology = mapping_row(d, offset)
            expected = expected_h(d, offset)
            ok = cohomology == expected
            all_ok &= ok
            print(
                f"{offset:2d} | {dims!s:>14} | {dminus.display():>14} | "
                f"{dzero.display():>14} | {ranks!s:>6} | {cohomology!s:>9} | "
                f"{expected!s:>9} | {'OK' if ok else 'FAIL'}"
            )
        print()
    if not all_ok:
        raise SystemExit("Characteristic-two verification failed.")
    print("ALL CHARACTERISTIC-TWO CHECKS PASSED for d=3,4,5.")


if __name__ == "__main__":
    main()
