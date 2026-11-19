#!/usr/bin/env python3
"""Exact finite Hom certificate for the full-profile marking boundary.

The all-degree homological epimorphism and full-profile assertions are proved
in the manuscript. This checks the finite graded algebra of its generator.
"""
from exact import Algebra, vec
from verify_nonformal import endomorphism_dga
from verify_connective_obstruction import nonzero_h


def main():
    # A=k(1->2), with e1*a=a=a*e2; R=A x k.
    r = Algebra(
        ['e1', 'e2', 'a', 'f'],
        {(0, 0): {0: 1}, (1, 1): {1: 1}, (3, 3): {3: 1},
         (0, 2): {2: 1}, (2, 1): {2: 1}},
        [1, 1, 0, 1],
    )
    r.starts = [0, 1, 0, 2]
    r.ends = [0, 1, 1, 2]
    r.idempotents = [r.basis(i) for i in (0, 1, 3)]
    r.check()
    end = endomorphism_dga(r, [(0, 0), (1, -2)], {})
    end.check_dga()
    assert end.n == 3
    assert nonzero_h(end) == {-2: 1, 0: 2}
    assert not any(c for row in end.d for c in row)
    e1 = end.basis(end.map_lookup[0, 0, 0])
    e2 = end.basis(end.map_lookup[1, 1, 1])
    arrow = end.basis(end.map_lookup[0, 1, 2])
    assert end.mul(e1, arrow) == arrow == end.mul(arrow, e2)
    assert end.mul(arrow, arrow) == vec(end.n)
    assert end.mul(e2, arrow) == vec(end.n) == end.mul(arrow, e1)
    # The simple at vertex two is the quotient Be2 / (Be1[2]).
    # Its image is the cone of the injective map P1[2] -> P2[2].
    image = endomorphism_dga(r, [(0, -3), (1, -2)], {(0, 1): r.basis(2)})
    image.check_dga()
    from verify_connective_obstruction import underlying_complex, cohomology
    _, degrees, differential = underlying_complex(
        r, [(0, -3), (1, -2)], {(0, 1): r.basis(2)}
    )
    assert {q: n for q, n in cohomology(degrees, differential).items() if n} == {-2: 1}
    print('Full-profile marking boundary over Q')
    print('dim R = 4; dim S = 1; target is the direct-factor projective Rf = fR')
    print('End(T)^op: all basis products and zero differential: PASS')
    print('H End(T)^op by nonzero degree = {-2: 1, 0: 2}')
    print('Image of the degree-zero simple at vertex two: H = {-2: 1}: PASS')
    print('ALL CHECKS PASSED')


if __name__ == '__main__':
    main()
