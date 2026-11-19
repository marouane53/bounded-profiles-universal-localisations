"""Small exact rational linear algebra. No third-party dependencies."""
from fractions import Fraction


def rref(rows, ncols=None):
    a = [[Fraction(x) for x in row] for row in rows]
    n = len(a[0]) if a else (ncols or 0)
    if any(len(row) != n for row in a):
        raise ValueError('Ragged matrix')
    pivots = []
    r = 0
    for c in range(n):
        p = next((i for i in range(r, len(a)) if a[i][c]), None)
        if p is None:
            continue
        a[r], a[p] = a[p], a[r]
        z = a[r][c]
        a[r] = [x / z for x in a[r]]
        for i in range(len(a)):
            if i != r and a[i][c]:
                z = a[i][c]
                a[i] = [x - z*y for x, y in zip(a[i], a[r])]
        pivots.append(c)
        r += 1
        if r == len(a):
            break
    return a, pivots


def rank(rows, ncols=None):
    return len(rref(rows, ncols)[1])


def nullspace(rows, ncols):
    a, p = rref(rows, ncols)
    out = []
    for f in range(ncols):
        if f in p:
            continue
        v = [Fraction(0)] * ncols
        v[f] = Fraction(1)
        for i, c in enumerate(p):
            v[c] = -a[i][f]
        out.append(v)
    return out


def column_rank(columns, nrows):
    return rank([[c[i] for c in columns] for i in range(nrows)], len(columns))
