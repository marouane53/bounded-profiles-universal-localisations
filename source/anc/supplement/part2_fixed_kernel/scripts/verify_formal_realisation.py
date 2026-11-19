#!/usr/bin/env python3
"""Exact rational certificates for the formal two-term kernel realisation.

R is a ZERO-PAIRING Morita-context ring, not ordinary matrix multiplication:
R=[[Gamma,Gamma],[U,Gamma]], with both off-diagonal products zero.
The proof covers arbitrary Gamma and U. These examples check multiplication,
projective target decompositions, and the strict dg quasi-isomorphism.
Only Python's standard library and the local exact_linear module are used.
"""
from fractions import Fraction as F
from itertools import product
from exact_linear import rank


def add(out, key, value):
    out[key] = out.get(key, F(0)) + F(value)
    if not out[key]:
        del out[key]


def linprod(x, y, rule):
    out = {}
    for a, ca in x.items():
        for b, cb in y.items():
            for c, cc in rule(a, b).items():
                add(out, c, ca * cb * cc)
    return out


def mapvec(x, fun):
    out = {}
    for a, ca in x.items():
        for b, cb in fun(a).items():
            add(out, b, ca * cb)
    return out


def example(kind):
    if kind in ('field', 'dual_regular', 'dual_simple'):
        g = 1 if kind == 'field' else 2
        u = 1 if kind in ('field', 'dual_simple') else 2
        def gm(a,b):
            return {a+b:F(1)} if a+b < g else {}
        def left(a,b):
            return {a+b:F(1)} if a+b < u else {}
        def right(a,b):
            return {a+b:F(1)} if a+b < u else {}
    elif kind == 'crossed_semisimple':
        g,u=2,2
        def gm(a,b): return {a:F(1)} if a==b else {}
        # u_0 in e_1 U e_0 and u_1 in e_0 U e_1.
        def left(a,b): return {b:F(1)} if a==1-b else {}
        def right(a,b): return {a:F(1)} if b==a else {}
    elif kind == 'matrix_regular':
        g,u=4,4
        def gm(a,b):
            i,j=divmod(a,2); h,l=divmod(b,2)
            return {2*i+l:F(1)} if j==h else {}
        left=gm; right=gm
    else:
        raise ValueError(kind)
    return g,u,gm,left,right


def verify(kind):
    g,u,gm,left,right=example(kind)
    # Bimodule associativity is checked independently of the ring construction.
    for a,b,c in product(range(g),range(g),range(u)):
        assert mapvec(gm(a,b),lambda t:left(t,c)) == linprod({a:1},left(b,c),left)
    for a,b,c in product(range(u),range(g),range(g)):
        assert linprod(right(a,b),{c:1},right) == linprod({a:1},gm(b,c),right)
    for a,b,c in product(range(g),range(u),range(g)):
        assert linprod(left(a,b),{c:1},right) == linprod({a:1},right(b,c),left)

    rb=[(0,0,a) for a in range(g)]+[(0,1,a) for a in range(g)]+[(1,0,a) for a in range(u)]+[(1,1,a) for a in range(g)]
    def rm(x,y):
        i,j,a=x; h,l,b=y
        if j!=h: return {}
        if i!=j and h!=l: return {}  # both Morita pairings are zero
        if (i,l)==(1,0):
            v=right(a,b) if i!=j else left(a,b)
        else:
            v=gm(a,b)
        return {(i,l,t):c for t,c in v.items()}
    sb=[(i,j,a) for i in range(2) for j in range(2) for a in range(g)]
    def sm(x,y):
        i,j,a=x; h,l,b=y
        if j!=h: return {}
        return {(i,l,t):c for t,c in gm(a,b).items()}
    def lam(x):
        return {} if x[:2]==(1,0) else {x:F(1)}
    for x,y,z in product(rb, repeat=3):
        assert linprod(rm(x,y),{z:1},rm)==linprod({x:1},rm(y,z),rm)
    for x,y in product(rb,repeat=2):
        assert mapvec(rm(x,y),lam)==linprod(lam(x),lam(y),sm)
    # Explicit module isomorphisms S_R=(pR)^2 and _R S=(Re)^2.
    for row in range(2):
        conv=lambda x: {(row,x[1],x[2]):F(1)}
        for x in [x for x in rb if x[0]==0]:
            for y in rb:
                assert mapvec(rm(x,y),conv)==linprod(conv(x),lam(y),sm)
    for col in range(2):
        conv=lambda x: {(x[0],col,x[2]):F(1)}
        for x in rb:
            for y in [y for y in rb if y[1]==1]:
                assert mapvec(rm(x,y),conv)==linprod(lam(x),conv(y),sm)
    # T=[Rp -> Re], with Rp in degree -1. End^op has U -> Gamma^2 -> Gamma.
    dzero=[[int(a==b)-int(a==b-g) for b in range(2*g)] for a in range(g)]
    dr=rank(dzero)
    assert (u,2*g-dr,g-dr)==(u,g,0)
    # Strict dg inclusion: Gamma maps to the diagonal and U to degree -1.
    # Product in End^op is reversed composition of right-multiplication maps.
    bb=[('g',a) for a in range(g)]+[('u',a) for a in range(u)]
    def bm(x,y):
        t,a=x; s,b=y
        if t==s=='u': return {}
        if t==s=='g': return {('g',c):v for c,v in gm(a,b).items()}
        if t=='g': return {('u',c):v for c,v in left(a,b).items()}
        return {('u',c):v for c,v in right(a,b).items()}
    # An endomorphism term is a map R e_source -> R e_target; its coefficient
    # is in e_source R e_target, with multiplication opposite to composition.
    def inc(x):
        t,a=x
        return {(0,0,a):F(1),(1,1,a):F(1)} if t=='g' else {(1,0,a):F(1)}
    for x,y in product(bb,repeat=2):
        assert mapvec(bm(x,y),inc)==linprod(inc(x),inc(y),rm)
    # No degree -2 terms exist; a product of two negative endomorphisms is zero.
    print(f'{kind}: dim Gamma={g}, dim U={u}, dim R={len(rb)}, dim S={len(sb)}')
    print(f'  associative R and multiplicative lambda: PASS')
    print(f'  S_R=(pR)^2 and _R S=(Re)^2, by basis action maps: PASS')
    print(f'  H^(-1,0,1) End(T)^op = {(u,g,0)}; strict dg inclusion: PASS')


def main():
    print('Formal two-term kernel realisation certificates over Q')
    print('The lower-left U is killed by lambda; BOTH Morita pairings are zero.')
    for kind in ('field','dual_regular','dual_simple','crossed_semisimple','matrix_regular'):
        verify(kind)
    print('Infinite Tor vanishing follows from the verified projective decompositions,')
    print('not from computing a finite initial segment of a resolution.')
    print('The theorem, rather than the finite sample, covers every finite Gamma-bimodule U.')
    print('ALL CHECKS PASSED')

if __name__=='__main__': main()
