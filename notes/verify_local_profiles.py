#!/usr/bin/env python3
"""Independent strict-bar, module-resolution, and incidence checks for the six-dimensional example."""

from collections import Counter
from itertools import product

from verify_profile_classification import PRIMES, clean, sparse_rank, strict_bar_differential, strict_bar_table


def all_words(n):
    if n < 0:
        return []
    return [''.join(w) for w in product('xz', repeat=n)]


def normal(word, from_right=False):
    while 'xxz' in word:
        i = word.rfind('xxz') if from_right else word.find('xxz')
        word = word[:i] + 'zxx' + word[i+3:]
    return word


def test_strict_bar():
    letters = {a: (0, 0, -1 if a in 'rw' else 0) for a in 'absrw'}
    strict = letters, {'aa': 's', 'br': 'w', 'rb': 'w'}, {'r': 's'}
    # Check the strict cochain signs before using any normal form.
    # u=s*, v=-r*, y=-w*, x=a*, z=b* in the convolution algebra.
    boundary = strict_bar_differential(*strict)

    def cup(left, right, word):
        if len(word) != 2:
            return 0
        first, second = word
        sign = (-1 if ((1-letters[first][2]) * (1-letters[second][2])) % 2 else 1)
        return sign * left.get(first, 0) * right.get(second, 0)

    x, z, v = {'a': 1}, {'b': 1}, {'r': -1}
    for length in (1, 2, 3):
        for chars in product(letters, repeat=length):
            word = ''.join(chars)
            bd = boundary(word)
            # d(u)=v-x²; d(y)=vz-zv. The dual signs are + for u and − for y.
            assert bd.get('s', 0) == v.get(word, 0) - cup(x, x, word)
            assert bd.get('w', 0) == cup(v, z, word) - cup(z, v, word)
    print('Strict bar transpose over the integers gives du=v−x² and dy=vz−zv.')
    observed = strict_bar_table(strict, 1, 6)
    expected = [sum('xxz' not in w for w in all_words(n)) for n in range(7)]
    assert [r[0, 0] for r in observed] == expected
    print('Six-dimensional strict dg algebra, direct bar dimensions through degree 6:', expected)
    formal = ({a: letters[a] for a in 'abw'}, {}, {})
    obs_formal = [r[0, 0] for r in strict_bar_table(formal, 1, 6)]
    assert obs_formal[:2] == [1, 2]
    assert all(obs_formal[n] == 2*obs_formal[n-1] + obs_formal[n-2] for n in range(2, 7))
    print('Formal comparison, direct bar dimensions through degree 6:', obs_formal)


def eval_word(word, even, odd, p):
    scalar = 1
    for n, a in enumerate(word):
        if a == 'z':
            scalar = scalar * (odd if n % 2 else even) % p
    return scalar


def test_resolutions():
    bases = [[w for w in all_words(n) if 'xxz' not in w] for n in range(9)]
    for n in range(11):
        for w in all_words(n):
            assert normal(w, True) == normal(w)
    for p in PRIMES:
        for even, odd in product(range(p), repeat=2):
            for n, basis in enumerate(bases):
                columns = []
                for step, prefix, scalar in [(1, 'z', even), (2, 'xz', odd)]:
                    if n < step:
                        continue
                    for suffix in bases[n-step]:
                        col = Counter({normal(prefix + suffix): 1})
                        col[normal('x' * step + suffix)] -= scalar
                        col = clean(col)
                        assert sum(v * eval_word(w, even, odd, p) for w, v in col.items()) % p == 0
                        columns.append(col)
                assert len(columns) == sparse_rank(columns, p) == len(basis) - 1
            # The degree-two commutator acts invertibly for unequal slopes.
            assert ((eval_word('zx', even, odd, p) - eval_word('xz', even, odd, p)) % p != 0) == (even != odd)
    print('Tested both rewrite orders through length 10 and all scalar/alternating slope resolutions through degree 8.')


def test_operator_detectors():
    # Two nilpotent operators with different products in the two orders.
    ops = (((0, 1), (0, 0)), ((0, 0), (1, 0)))
    bases = [[w for w in all_words(n) if 'xxz' not in w] for n in range(8)]
    for p in PRIMES:
        for n, basis in enumerate(bases):
            columns = []
            for step, prefix, matrix in [(1, 'z', ops[0]), (2, 'xz', ops[1])]:
                if n < step:
                    continue
                for suffix in bases[n-step]:
                    for v in range(2):
                        col = Counter({(v, normal(prefix + suffix)): 1})
                        for j in range(2):
                            col[j, normal('x'*step + suffix)] -= matrix[j][v]
                        columns.append(clean(col))
            assert len(columns) == sparse_rank(columns, p) == 2*len(basis)-2
            for col in columns:
                image = [0, 0]
                for (v, w), coeff in col.items():
                    vector = [int(j == v) for j in range(2)]
                    for j, a in enumerate(w):
                        if a == 'z':
                            op = ops[j % 2]
                            vector = [sum(a*b for a, b in zip(row, vector)) for row in op]
                    for j in range(2):
                        image[j] += coeff*vector[j]
                assert all(value % p == 0 for value in image)
    print('Checked a noncommuting two-by-two operator pair and its free resolution through degree 7.')


def divisible(f, divisor, p):
    remainder = list(f)
    while remainder and len(remainder) >= len(divisor):
        offset = len(remainder) - len(divisor)
        coeff = remainder[-1]
        for j, a in enumerate(divisor):
            remainder[offset+j] = (remainder[offset+j] - coeff*a) % p
        while remainder and remainder[-1] == 0:
            remainder.pop()
    return not remainder


def irreducibles(p):
    polys = []
    for degree in range(1, 4):
        for coeffs in product(range(p), repeat=degree):
            f = coeffs + (1,)
            if not any(divisible(f, g, p) for g in polys if 2*(len(g)-1) <= degree):
                polys.append(f)
    return polys


def multiply_T(f, i):
    d = len(f)-1
    return {i+1: 1} if i+1 < d else {j: -f[j] for j in range(d) if f[j]}


def test_closed_points():
    for p in PRIMES:
        polys = irreducibles(p)
        assert len(polys) == p + (p*p-p)//2 + (p**3-p)//3
        for f, g in product(polys, repeat=2):
            d, e = len(f)-1, len(g)-1
            # Multiplication by T⊗1 - 1⊗T on k_f ⊗ k_g, over the base field.
            # Transposing either companion matrix gives a similar matrix.
            columns = []
            for i, j in product(range(d), range(e)):
                col = Counter({(a, j): v for a, v in multiply_T(f, i).items()})
                for b, val in multiply_T(g, j).items():
                    col[i, b] -= val
                columns.append(clean(col))
            assert d*e - sparse_rank(columns, p) == (d if f == g else 0)
        print(f'F{p}: checked {len(polys)**2} ordered closed-point pairs, all monic irreducibles of degrees 1 to 3.')


if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('Enable Python assertions to run the verification.')
    test_strict_bar()
    test_resolutions()
    test_closed_points()
    test_operator_detectors()
    print('All finite checks passed. The note gives the proofs in arbitrary degree and characteristic.')
