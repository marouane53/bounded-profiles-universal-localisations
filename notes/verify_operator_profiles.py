#!/usr/bin/env python3
"""Independent finite checks for the profile comparison and its proof examples."""

from itertools import product

from compare_operator_pairs import Basis, Field, compare, profile_included, read_pair


def matrix_sum(a, b, sign=1):
    return [[x+sign*y for x, y in zip(r, t)] for r, t in zip(a, b)]


def multiply(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def identity(n, scale=1):
    return [[scale*int(i == j) for j in range(n)] for i in range(n)]


def rank(a, field):
    space = Basis(field)
    for row in a:
        space.add([field.reduce(x) for x in row])
    return len(space.rows)


def over_f2(a):
    return tuple(tuple(x % 2 for x in row) for row in a)


def independent_types_f2(pair, bases):
    # Enumerate actual invariant lines; in this dimension a non-simple pair
    # has two scalar factors. For a simple pair enumerate its whole GL_2 orbit.
    conjugates = []
    for basis, inverse in bases:
        conjugated = tuple(over_f2(multiply(inverse, multiply(op, basis))) for op in pair)
        if all(m[1][0] == 0 for m in conjugated):
            return frozenset(('scalar', tuple(sorted(m[i][i] for m in conjugated))) for i in range(2))
        conjugates.extend([conjugated, conjugated[::-1]])
    return frozenset([('simple', min(conjugates))])


def check_pair_comparison():
    field = Field(2)
    mats = [[[a,b], [c,d]] for a,b,c,d in product(range(2), repeat=4)]
    bases = [(m, [[m[1][1], m[0][1]], [m[1][0], m[0][0]]]) for m in mats
             if (m[0][0]*m[1][1]-m[0][1]*m[1][0]) % 2]
    assert len(bases) == 6
    cases = []
    representatives = {}
    for pair in product(mats, repeat=2):
        support = independent_types_f2(pair, bases)
        data = read_pair(list(pair), field)
        cases.append((support, data))
        representatives.setdefault(support, data)
    calls = 0
    for support, data in cases:
        for other, exemplar in representatives.items():
            assert profile_included(data, exemplar, field) == (other <= support)
            assert profile_included(exemplar, data, field) == (support <= other)
            calls += 2
    print(f'F2: independently classified all 256 pairs of two-by-two matrices into {len(representatives)} profiles.')
    print(f'Checked {calls} directed inclusions against invariant lines and changes of basis.')

    double = [[[0, 1], [0, 0]], [[0, 1], [0, 0]]]
    zero = [[[0]], [[0]]]
    assert compare({'first': zero, 'second': double})['equal']
    for prime in (2, 3, 5):
        empty = compare({'prime':prime, 'first':[[], []], 'second':zero})
        assert not empty['first_subset_of_second'] and empty['second_subset_of_first']
    # A nonrational simple over Q has its conjugate eigenvalues in one factor.
    rotation = [[0,-1],[1,0]]
    rational_simple = [rotation, rotation]
    assert compare({'first':rational_simple,'second':rational_simple})['equal']
    assert not any(compare({'first':rational_simple,'second':zero}).values())
    assert compare({'first':[[['1/2']],[[2]]], 'second':[[[2]],[['1/2']]]})['equal']
    print('Checked nonsplit extensions, empty supports, rational fractions and a nonsplit-field simple.')


def strip(word, letter):
    if not word or word[-1] != letter:
        return 0, ''
    return (-1 if (len(word)-1) % 2 else 1), word[:-1]


def fibre_check(limit=12):
    words = [['']]
    for n in range(1, limit+1):
        words.append([c+w for w in words[-1] for c in 'ab' if not (c == 'a' and w.startswith('a'))])
    def basis_r(n):
        return [] if n < 0 else ([''] if n == 0 else (['a','b'] if n == 1 else ['a'*n]))
    def times(r, c):
        return r+c if not r or (c == 'a' and 'b' not in r) else None
    for weight in range(limit+1):
        dimensions, ranks = [], []
        for n in range(weight+1):
            bases = [(r,w) for r in basis_r(weight-n) for w in words[n]]
            dimensions.append(len(bases))
            images = set()
            for r,w in bases:
                if not w:
                    continue
                value = times(r, w[0])
                if value is not None:
                    images.add((value, w[1:]))
            ranks.append(len(images))  # Each nonzero column is a single basis vector.
        ranks.append(0)
        assert [d-ranks[i]-ranks[i+1] for i,d in enumerate(dimensions)] == [int(weight == 0)]+[0]*weight
    for group in words[1:]:
        for w in group:
            for letter in 'ab':
                before, tail = strip(w, letter)
                after, suffix = strip(w[1:], letter)
                if before and tail:
                    assert after == -before and tail[1:] == suffix
                else:
                    assert not after
    print(f'Special fibre: verified exactness in internal weights 0 to {limit} and the signed stripping cocycles.')


def concrete_separators():
    u = [[2,0,0],[0,3,0],[0,0,4]]
    v = [[0,0,1],[1,0,0],[0,1,0]]
    transposed = [list(row) for row in zip(*v)]
    def separator(a, b):
        comm = matrix_sum(multiply(a,b),multiply(b,a),-1)
        cube = multiply(comm,multiply(comm,comm))
        q = identity(3)
        for t in (2,3,4):
            q = multiply(q,matrix_sum(a,identity(3,t),-1))
        return matrix_sum(matrix_sum(cube,identity(3,2),-1),q)
    field = Field()
    assert rank(separator(u,v), field) == rank(separator(v,u), field) == 3
    assert rank(separator(u,transposed),field) == 0
    comparison = compare({'first':[u,v], 'second':[u,transposed]})
    assert not any(comparison.values())
    print('Over Q the transpose example has identical scalar-pencil determinants and a separating degree-six polynomial.')

    j, k = [[0,1],[0,0]], [[0,0],[1,0]]
    for prime in (2,3,5):
        assert rank(matrix_sum(multiply(j,k),multiply(k,j),-1),Field(prime)) == 2
    print('The commutator has empty scalar domain and an invertible two-dimensional simple evaluation in all three test fields.')


def check_schur_witnesses():
    # Test the finite-word graph used to compare normalized matrix pencils.
    polys = [{'UV':1}, {'UV':1,'VU':-1}, {'U':1,'VUV':1,'UU':-2}]
    for p in (2,3,5):
        field = Field(p)
        candidates = [[[0,1],[0,0]], [[0,0],[1,0]], [[1,1],[0,1]], [[0,1],[-1,0]]]
        for poly in polys:
            vertices = 1
            edges = []
            for word, coeff in poly.items():
                visit = [0] + list(range(vertices,vertices+len(word)-1)) + [0]
                vertices += len(word)-1
                edges.extend((visit[i],visit[i+1],letter,coeff if i == 0 else 1)
                             for i,letter in enumerate(word))
            for a,b in product(candidates, repeat=2):
                ops = {'U':a,'V':b}
                linear = identity(2*vertices)
                result = identity(2)
                for i,j,c,n in edges:
                    for s,t in product(range(2),repeat=2):
                        linear[2*i+s][2*j+t] -= n*ops[c][s][t]
                for word,n in poly.items():
                    term = identity(2)
                    for c in word:
                        term = multiply(term,ops[c])
                    result = matrix_sum(result,term,-n)
                assert 2*vertices-rank(linear,field) == 2-rank(result,field)
    print('The path Schur complements preserve nullity for 144 polynomial/tuple/field combinations.')


if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('Enable Python assertions to run the checks.')
    check_pair_comparison()
    fibre_check()
    concrete_separators()
    check_schur_witnesses()
    print('All checks passed. The note contains the unrestricted proofs.')
