#!/usr/bin/env python3
"""Exact finite identities underlying the coordinate-symmetry example.

Uses integer basis multiplication and ranks over Q and F_2. The manuscript
proves the all-degree statements from the displayed periodic resolution;
this script checks its repeating step, not an arbitrary finite Ext cutoff.
"""
from fractions import Fraction

BASIS = ((0, 0), (1, 0), (0, 1), (1, 1))  # 1, x, y, xy
SWAP = (0, 2, 1, 3)


def product(i, j):
    a, b = BASIS[i][0] + BASIS[j][0], BASIS[i][1] + BASIS[j][1]
    return BASIS.index((a, b)) if a < 2 and b < 2 else None


def matrix_for(variable):
    return [[int(product(variable, col) == row) for col in range(4)]
            for row in range(4)]


def multiply(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def rank(a, prime=None):
    a = [[Fraction(x) if prime is None else x % prime for x in row] for row in a]
    pivot_row = 0
    for col in range(len(a[0])):
        row = next((i for i in range(pivot_row, len(a)) if a[i][col]), None)
        if row is None:
            continue
        a[pivot_row], a[row] = a[row], a[pivot_row]
        inverse = 1 / a[pivot_row][col] if prime is None else pow(a[pivot_row][col], -1, prime)
        a[pivot_row] = [x * inverse if prime is None else x * inverse % prime
                        for x in a[pivot_row]]
        for i in range(len(a)):
            if i == pivot_row:
                continue
            c = a[i][col]
            a[i] = [x - c * y if prime is None else (x - c * y) % prime
                    for x, y in zip(a[i], a[pivot_row])]
        pivot_row += 1
        if pivot_row == len(a):
            break
    return pivot_row


def main():
    assert SWAP[0] == 0 and all(SWAP[SWAP[i]] == i for i in range(4))
    for i in range(4):
        for j in range(4):
            ij = product(i, j)
            swapped = None if ij is None else SWAP[ij]
            assert product(SWAP[i], SWAP[j]) == swapped
    x, y = matrix_for(1), matrix_for(2)
    zero4 = [[0] * 4 for _ in range(4)]
    assert multiply(x, x) == multiply(y, y) == zero4
    zero2, jordan = [[0, 0], [0, 0]], [[0, 0], [1, 0]]
    # On M=C/(x): x=0, y=J. On N=C/(y): x=J, y=0.
    action_m, action_n = (zero2, jordan), (jordan, zero2)
    assert action_m[::-1] == action_n and multiply(jordan, jordan) == zero2
    for prime in (None, 2):
        assert rank(x, prime) == rank(y, prime) == 2
        # Square zero plus half rank gives image=kernel, hence the repeat is exact.
        assert rank(jordan, prime) == 1
        assert 2 - 2 * rank(zero2, prime) == 2  # self-Ext in positive degrees
        assert 2 - 2 * rank(jordan, prime) == 0  # cross-Ext in positive degrees
    print('Coordinate symmetry for C = k[x,y]/(x^2,y^2)')
    print('Swap x <-> y: unital involution; all 16 multiplication pairs: PASS')
    print('Module pullback exchanges C/(x) and C/(y): PASS')
    print('Repeating free-resolution step: image = kernel over Q and F_2: PASS')
    print('Repeating Hom step: self-Ext dimension 2; cross-Ext dimension 0: PASS')
    print('All-degree Ext statements use the periodicity proof in the manuscript.')
    print('ALL COORDINATE-SYMMETRY CHECKS PASSED')


if __name__ == '__main__':
    main()
