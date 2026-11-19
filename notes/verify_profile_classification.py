#!/usr/bin/env python3
"""Exact finite checks for the profile notes. Uses only the Python standard library."""

from collections import Counter, defaultdict
from functools import cache
from itertools import product

PRIMES = (2, 3, 5)


def add_multiple(total, column, coefficient):
    for key, value in column.items():
        total[key] += coefficient * value


def clean(column):
    return {key: value for key, value in column.items() if value}


class ColumnBasis:
    """Sparse column elimination in a prime field, using word labels as row indices."""

    def __init__(self, prime):
        self.prime = prime
        self.rows = {}

    def insert(self, column):
        p = self.prime
        v = {i: a % p for i, a in column.items() if a % p}
        while v:
            i = max(v)
            a = v[i]
            if i not in self.rows:
                inv = pow(a, -1, p)
                self.rows[i] = {j: b * inv % p for j, b in v.items()}
                return True
            for j, b in self.rows[i].items():
                value = (v.get(j, 0) - a * b) % p
                if value:
                    v[j] = value
                else:
                    v.pop(j, None)
        return False


def sparse_rank(columns, prime):
    basis = ColumnBasis(prime)
    for column in columns:
        basis.insert(column)
    return len(basis.rows)


def local_table():
    return {'t': (0, 0, 0), 'c': (0, 0, 0), 'r': (0, 0, -1),
            'w': (0, 0, -1)}, {'tt': 'c', 'tr': 'w'}, {'r': 'c'}


def triangle_table(with_bridge):
    letters = {'a': (0, 1, 0), 'b': (1, 2, 0), 'c': (2, 0, 0),
               's': (1, 0, 0), 'r': (1, 0, -1), 'w': (0, 0, -1)}
    if with_bridge:
        letters.update({'h': (0, 3, 0), 'l': (3, 3, 0)})
    return letters, {'bc': 's', 'ar': 'w'}, {'r': 's'}


def verify_table(letters, products, differential):
    def mul(a, b):
        return products.get(a + b) if a and b else None

    def singleton(a, sign=1):
        return Counter({a: sign}) if a else Counter()

    for a, b in products:
        src, _, deg = letters[a]
        _, end, other_deg = letters[b]
        assert letters[a][1] == letters[b][0]
        assert letters[mul(a, b)] == (src, end, deg + other_deg)
    for a, da in differential.items():
        src, end, deg = letters[a]
        assert letters[da] == (src, end, deg + 1) and da not in differential
    for a, b, c in product(letters, repeat=3):
        assert mul(mul(a, b), c) == mul(a, mul(b, c))
    for a, b in product(letters, repeat=2):
        rhs = singleton(mul(differential.get(a), b))
        rhs.update(singleton(mul(a, differential.get(b)), (-1 if letters[a][2] % 2 else 1)))
        assert clean(rhs) == clean(singleton(differential.get(mul(a, b))))


def strict_bar_differential(letters, products, differential):
    @cache
    def boundary(word):
        answer = Counter()
        sign = 1
        for i, a in enumerate(word):
            deg = letters[a][2]
            if a in differential:
                answer[word[:i] + differential[a] + word[i+1:]] -= sign
            pair = word[i:i+2]
            if pair in products:
                answer[word[:i] + products[pair] + word[i+2:]] += sign * (-1 if deg % 2 else 1)
            if (1 - deg) % 2:
                sign = -sign
        return clean(answer)

    return boundary


def strict_bar_table(table, vertices, limit):
    """Homology dimensions by both endpoints of the original strictly associative dg algebra."""
    letters, products, differential = table
    verify_table(*table)

    boundary = strict_bar_differential(letters, products, differential)

    # Degree n in the dual cochains corresponds to bar words of weight n.
    paths = [defaultdict(list) for _ in range(limit + 2)]
    for i in range(vertices):
        paths[0][i, i].append('')
    ranks = {p: [] for p in PRIMES}
    for n in range(limit + 2):
        for p in PRIMES:
            ranks[p].append({})
        for (start, end), words in paths[n].items():
            columns = [boundary(w) for w in words]
            for col in columns:
                square = Counter()
                for w, val in col.items():
                    add_multiple(square, boundary(w), val)
                assert not clean(square)
            for p in PRIMES:
                ranks[p][n][start, end] = sparse_rank(columns, p)
            for a, (src, dst, deg) in letters.items():
                next_n = n + 1 - deg
                if end == src and next_n < len(paths):
                    paths[next_n][start, dst].extend(w + a for w in words)
    dimensions = {}
    for p in PRIMES:
        dimensions[p] = [
            {(i, j): len(paths[n].get((i, j), ()))
             - ranks[p][n].get((i, j), 0) - ranks[p][n+1].get((i, j), 0)
             for i, j in product(range(vertices), repeat=2)}
            for n in range(limit + 1)]
    assert all(dimensions[p] == dimensions[2] for p in PRIMES)
    return dimensions[2]


def local_numbers(limit):
    result = []
    for n in range(limit + 1):
        value = int(n <= 2)
        for shift in (3, 4):
            if n >= shift:
                value += result[n-shift]
        result.append(value)
    return result


def check_strict_models():
    local = strict_bar_table(local_table(), 1, 7)
    assert [row[0, 0] for row in local] == local_numbers(7)
    print('Five-dimensional strict dg model, direct bar ranks through degree 7:', local_numbers(7))

    for vertices in (3, 4):
        table = triangle_table(vertices == 4)
        observed = strict_bar_table(table, vertices, 9)
        # For left simples, paths from j to i compute Ext(S_i,S_j).
        resolutions = [[0, 2, 1], [1, 0, 2, 1], [2, 1, 0, 2, 1]]
        for n, row in enumerate(observed):
            for i, seq in enumerate(resolutions):
                for j in range(vertices):
                    assert row[j, i] == int(n < len(seq) and seq[n] == j)
            if vertices == 4:
                # The radical of P_3 is the direct sum S_0 and S_3.
                for j in range(vertices):
                    expected = (int(j == 3) if n == 0 else
                                observed[n-1][j, 0] + observed[n-1][j, 3])
                    assert row[j, 3] == expected
        totals = [sum(row.values()) for row in observed]
        print(f'{vertices}-vertex nonformal model, direct bar totals through degree 9:', totals)

        letters, _, _ = table
        formal = ({a: v for a, v in letters.items() if a not in ('r', 's')}, {}, {})
        f = strict_bar_table(formal, vertices, 9)
        assert all(row[0, 0] > 0 for n, row in enumerate(f) if n % 3 == 0)
        if vertices == 4:
            assert all(row[3, 3] == 1 and all(row[3, i] == 0 for i in range(3)) for row in f)
        print(f'{vertices}-vertex formal model, direct bar totals through degree 9:',
              [sum(row.values()) for row in f])


@cache
def words(n):
    if n < 0:
        return ()
    if n == 0:
        return ('',)
    return tuple('x' + w for w in words(n-1)) + tuple('y' + w for w in words(n-2))


@cache
def deriv(word):
    coefficients = Counter()
    sign = 1
    for i, letter in enumerate(word):
        if letter == 'y':
            coefficients[word[:i] + 'xxx' + word[i+1:]] += sign
        else:
            sign = -sign
    return clean(coefficients)


@cache
def normal_cycles(n):
    result = [{'x' * n: 1}] if n in (0, 1, 2) else []
    for step, cycle in [(3, {'xy': 1, 'yx': 1}), (4, {'xxy': 1, 'yxx': -1})]:
        if n < step:
            continue
        for suffix in normal_cycles(n-step):
            coefficients = Counter()
            for a, x in cycle.items():
                for b, y in suffix.items():
                    coefficients[a + b] += x * y
            result.append(clean(coefficients))
    return result


def check_normal_forms(limit=16):
    expected = local_numbers(limit)
    for p in PRIMES:
        for n in range(limit + 1):
            cycles = normal_cycles(n)
            image = ColumnBasis(p)
            for word in words(n-1):
                image.insert(deriv(word))
            kernel_dimension = len(words(n)) - sparse_rank([deriv(w) for w in words(n)], p)
            assert kernel_dimension - len(image.rows) == expected[n] == len(cycles)
            for cycle in cycles:
                delta = Counter()
                for w, coefficient in cycle.items():
                    add_multiple(delta, deriv(w), coefficient)
                assert not clean(delta) and image.insert(cycle)
            assert len(image.rows) == kernel_dimension
    print(f'Bar-dual normal forms are a cohomology basis through degree {limit}:', expected)


if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('Enable Python assertions to run the verification.')
    check_strict_models()
    check_normal_forms()
    print('All degreewise rank and basis checks passed over F2, F3, and F5.')
