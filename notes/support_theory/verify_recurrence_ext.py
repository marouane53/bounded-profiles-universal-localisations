"""Compute both sides of rank-one resolutions and pair Ext from the multiplication table."""
from itertools import product
from math import comb

from check_lib import rank, require_asserts


def apply_q(q, v, p):
    return tuple(sum(x*y for x, y in zip(row, v)) % p for row in q)


def products_of_basis(q, p):
    r = len(q)
    size = 2*r + 2
    # 0 is the unit, 1 is a, indices 2..r+1 are the b_i, the rest the w_i.
    table = {(0, i): {i: 1} for i in range(size)}
    table.update({(i, 0): {i: 1} for i in range(size)})
    for i in range(r):
        table[1, 2+i] = {2+r+i: 1}
        table[2+i, 1] = {2+r+j: -q[j][i] % p for j in range(r)}
    return table, size


def product_maps(q, v, table, size, p):
    l = (0, 1, *v, *(0 for _ in q))
    left, right = [], []
    for i in range(size):
        col_l, col_r = [0]*size, [0]*size
        for j, coeff in enumerate(l):
            if not coeff:
                continue
            for out, value in table.get((j, i), {}).items():
                col_l[out] += coeff*value
            for out, value in table.get((i, j), {}).items():
                col_r[out] += coeff*value
        left.append(tuple(x % p for x in col_l))
        right.append(tuple(x % p for x in col_r))
    return left, right


def check_composition(first, second, p):
    for column in first:
        result = [0] * len(column)
        for coeff, image in zip(column, second):
            if coeff:
                for i, value in enumerate(image):
                    result[i] += coeff*value
        assert all(x % p == 0 for x in result)


def all_gl(r, p):
    for entries in product(range(p), repeat=r*r):
        q = [entries[i*r:(i+1)*r] for i in range(r)]
        if rank(q, p) == r:
            yield q


def verify_matrix(q, p, steps=9):
    r = len(q)
    table, size = products_of_basis(q, p)
    vectors = list(product(range(p), repeat=r))
    maps = {v: product_maps(q, v, table, size, p) for v in vectors}
    for v, (left, right) in maps.items():
        qv = apply_q(q, v, p)
        next_left, next_right = maps[qv]
        assert rank(left, p) == rank(right, p) == r+1
        check_composition(next_left, left, p)
        check_composition(right, next_right, p)
    cases = 0
    for v in vectors:
        # Compute the entries of the differentials in the ordinary algebra
        # from the first columns of the successive regular-module maps.
        current = v
        coefficients = []
        for _ in range(steps+1):
            coefficients.append(maps[current][0][0])
            current = apply_q(q, current, p)
        for ell in vectors:
            for c in range(p):
                # Independent computation: evolve the functional by right
                # multiplication by Q, preserving the original vector v.
                row = ell
                indicators = []
                ranks = []
                for column in coefficients:
                    value = (c + sum(i*j for i, j in zip(row, v))) % p
                    indicators.append(int(value == 0))
                    row = tuple(sum(row[i]*q[i][j] for i in range(r)) % p
                                for j in range(r))
                    action = (c*column[1] + sum(x*y for x, y in zip(ell, column[2:2+r]))) % p
                    # The two-dimensional module has e*l=action*f and f*l=0.
                    ranks.append(int(bool(action)))
                for j in range(1, steps+1):
                    assert 2 - ranks[j-1] - ranks[j] == indicators[j-1] + indicators[j]
                cases += 1
    return cases


def sparse_example():
    for p in (2, 3, 5, 7):
        # Compare the ordinary example's scalar to the basis coefficients
        # of (1+t)^j-t^j-1, directly over the polynomial ring.
        zeros = [False] + [all(comb(j, i) % p == 0 for i in range(1, j)) for j in range(1, 401)]
        powers = set()
        j = 1
        while j <= 400:
            powers.add(j)
            j *= p
        assert {j for j in range(401) if zeros[j]} == powers
        dims = [1] + [int(zeros[j-1]) + int(zeros[j]) for j in range(1, 401)]
        assert dims == [1] + [int(j-1 in powers) + int(j in powers) for j in range(1, 401)]


if __name__ == "__main__":
    require_asserts()
    cases = matrices = 0
    for r, p in ((1, 2), (1, 3), (2, 2), (2, 3), (3, 2)):
        for q in all_gl(r, p):
            cases += verify_matrix(q, p)
            matrices += 1
    print(f"{matrices} invertible Q, {cases} ordinary module pairs, Ext degrees 1–9.")
    print("Checked ranks and compositions of the regular free complex and its algebra dual.")
    sparse_example()
    print("Sparse ordinary-module example through degree 400 over F_p(t), p=2,3,5,7.")
