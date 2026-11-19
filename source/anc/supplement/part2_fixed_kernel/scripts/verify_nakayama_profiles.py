#!/usr/bin/env python3
"""Exact rational stable-Hom and profile certificates for d=3,4.

Right A_d modules M(i,l) have factors S_i,...,S_(i+l-1).
Hom is computed from arrow-commutativity equations; maps through all
indecomposable projectives are quotiented out. Symbolic formulas and the
proof in the Phase 1 text cover every shift; this checks a full period.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from exact_linear import nullspace, column_rank


@lru_cache(None)
def hom(d, x, y):
    n = 2*d
    i, l = x
    j, m = y
    sx = {(i+r) % n for r in range(l)}
    sy = {(j+r) % n for r in range(m)}
    common = sorted(sx & sy)
    ix = {v: a for a, v in enumerate(common)}
    ax = {(i+r) % n for r in range(l-1)}
    ay = {(j+r) % n for r in range(m-1)}
    equations = []
    for v in range(n):
        row = [0]*len(common)
        if v in ay and v in ix:
            row[ix[v]] += 1
        vp = (v+1) % n
        if v in ax and vp in ix:
            row[ix[vp]] -= 1
        equations.append(row)
    ans = []
    for z in nullspace(equations, len(common)):
        vec = [F(0)] * n
        for v, a in ix.items():
            vec[v] = z[a]
        ans.append(tuple(vec))
    return tuple(ans)


@lru_cache(None)
def stable_data(d, x, y):
    h = hom(d, x, y)
    factors = []
    for i in range(2*d):
        p = (i, d)
        for f in hom(d, x, p):
            for g in hom(d, p, y):
                factors.append(tuple(a*b for a,b in zip(f,g)))
    r = column_rank(factors, 2*d)
    return len(h)-r, len(h), r


def omega(d, x):
    i,l = x
    return ((i+l) % (2*d), d-l)


def suspension(d, x):
    i,l = x
    return ((i+l-d) % (2*d), d-l)


def orbit(d, x):
    ans = []
    while x not in ans:
        ans.append(x)
        x = suspension(d,x)
    assert len(ans)==4
    return tuple(ans)


def total(d,x,y):
    return tuple(stable_data(d,x,z)[0] for z in orbit(d,y))


def name(x):
    i,l=x
    return f'S{i+1}' if l==1 else f'M({i+1},{l})'


def main():
    print('Nakayama stable-Hom/profile certificate over Q')
    for d in (3,4):
        modules=[(i,l) for i in range(2*d) for l in range(1,d)]
        for x in modules:
            for y in modules:
                r=(x[0]-y[0])%(2*d)
                formula=int(r<y[1] and y[1]<=r+x[1]<d)
                assert stable_data(d,x,y)[0]==formula
        print(f'd={d}: all {len(modules)**2} stable-Hom entries match r<m<=r+l<d: PASS')
        for x in modules:
            z=x
            for _ in range(4): z=omega(d,z)
            assert z==x
        print(f'd={d}: Omega^4 fixes all indecomposables: PASS')
        for j in range(d):
            o=orbit(d,(j,1))
            mat=[[stable_data(d,x,y)[0] for y in o] for x in o]
            assert mat==[[int(a==b) for b in range(4)] for a in range(4)]
        print(f'd={d}: all singleton suspension orbits have identity Hom^0 matrix: PASS')
    d=4
    reps=[(0,1),(1,1),(0,2)]
    print('\nd=4, J={1,2}; row=source, column=detector; entries are dimensions in shifts 0,1,2,3')
    for x in reps:
        print(name(x), [total(d,x,y) for y in reps])
    expected=[[1,1,0],[0,1,1],[1,0,1]]
    actual=[[int(any(total(d,x,y))) for y in reps] for x in reps]
    assert actual==expected,(actual,expected)
    for a,x in enumerate(reps):
        for b,y in enumerate(reps):
            for r,xx in enumerate(orbit(d,x)):
                for s,yy in enumerate(orbit(d,y)):
                    assert stable_data(d,xx,yy)[0]==total(d,x,y)[(s-r)%4]
    print('All 144 degree-zero Hom entries recovered from the representative shift table: PASS')
    print('Non-split triangle: S2 -> M(1,2) -> S1 -> S2[1]')
    print('Realising W for [full, orbit(S1), orbit(S2), orbit(M(1,2)), zero]:')
    print('  [0, D M(1,2), D S1, D S2, D(S1+S2)]')
    for j in range(4):
        reps=[(j,1),((j+1)%8,1),(j,2)]
        expected_members=[{0,1,2},{0},{1},{2},set()]
        detectors=[[],[reps[2]],[reps[0]],[reps[1]],[reps[0],reps[1]]]
        for ds,want in zip(detectors,expected_members):
            got={a for a,x in enumerate(reps) if not any(any(total(4,x,y)) for y in ds)}
            assert got==want,(j,ds,got,want)
    print('All 4 cyclic adjacent J and all 5 profiles: PASS')
    for j in range(2):
        reps=[(j,1),((j+2)%8,1)]
        for a in range(2):
            for b in range(2):
                assert bool(any(total(4,reps[a],reps[b])))==(a==b)
        for keep in (set(),{0},{1},{0,1}):
            ds=[reps[a] for a in range(2) if a not in keep]
            got={a for a,x in enumerate(reps) if not any(any(total(4,x,y)) for y in ds)}
            assert got==keep
    print('All 2 opposite J and all 4 profiles: PASS')
    print('\nFull list of nonempty admissible J (1-based labels):')
    for d in (3,4):
        for size in range(1,d-1):
            for J in combinations(range(1,d+1),size):
                if size==1: kind,count='singleton',2
                elif (J[1]-J[0])%d in (1,d-1): kind,count='adjacent',5
                else: kind,count='opposite',4
                print(f'  d={d}, J={J}: {kind}, {count} thick subcategories, all realisable')
    print('ALL CHECKS PASSED')

if __name__=='__main__': main()
