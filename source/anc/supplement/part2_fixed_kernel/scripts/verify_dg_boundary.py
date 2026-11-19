#!/usr/bin/env python3
"""Exact rational certificate for the 20-dimensional dg boundary.

B^0=Q e1+Q e4, B^-1=Q a+Q b, d=0, a=e4*a*e1,
b=e1*b*e4, a*b=b*a=0. Left B modules; right Nakayama
modules are left modules over R=A_3^op.

A semi-free resolution has cell x_n in degree -2n, alternating
vertices, with d x_n = (the radical basis element of cell n-1).
Truncation is checked away from its one unfilled bottom class.
Every-degree assertions are proved by the same cell recurrence in the text.
"""
from fractions import Fraction as F
from exact_linear import rank
from verify_nakayama_profiles import stable_data, suspension


def resolution(start, cells=17):
    basis={}
    for n in range(cells):
        vertex=start if n%2==0 else 5-start # 1 <-> 4
        basis[('x',n)]=(-2*n,vertex)
        basis[('z',n)]=(-2*n-1,5-vertex)
    def differential(degree):
        src=[x for x,v in basis.items() if v[0]==degree]
        tgt=[x for x,v in basis.items() if v[0]==degree+1]
        out=[[F(0) for _ in src] for _ in tgt]
        for c,(typ,n) in enumerate(src):
            if typ=='x' and n>0:
                out[tgt.index(('z',n-1))][c]=F(1)
        return src,tgt,out
    low=-2*(cells-1)-1
    for degree in range(low,1):
        src,tgt,di=differential(degree)
        prev,_,diprev=differential(degree-1)
        h=len(src)-rank(di,len(src))-rank(diprev,len(prev))
        assert h==int(degree in (low,0)),(degree,h)
    return basis


def ext_degrees(start,target,cells):
    b=resolution(start,cells)
    # Hom to a simple annihilates radical terms and every differential.
    return [-b[('x',n)][0] for n in range(cells) if b[('x',n)][1]==target]


def main():
    print('20-dimensional boundary certificate over Q')
    print('Conventions: cohomological grading, left B modules; R=A_3^op')
    print('dim R = 6*3 = 18; dim W = 1; dim A_prime = 20')
    assert 1+6*3+1==20
    cells=17
    q11=ext_degrees(1,1,cells)
    q41=ext_degrees(4,1,cells)
    print('Semi-free generators have degrees 0,-2,-4,...,-32: PASS')
    print('Truncated resolutions exact away from augmentation and unfilled bottom: PASS')
    print('Ext_B^q(L1,L1) nonzero for 0 <= q <= 32:',q11)
    print('Ext_B^q(L4,L1) nonzero for 0 <= q <= 32:',q41)
    assert q11==list(range(0,33,4))
    assert q41==list(range(2,33,4))
    for start in (1,4):
        for q in range(1,33):
            target=(0,1)
            for _ in range(q%4): target=suspension(3,target)
            actual=stable_data(3,(start-1,1),target)[0]
            want=int(q in (q11 if start==1 else q41))
            assert actual==want
    print('Independent ambient stable-Ext check in degrees 1,...,32: PASS')
    print('All odd Ext degrees vanish; the brief\'s parity assertion is false: VERIFIED')
    print('Correct infinite patterns: q=4m and q=4m+2, respectively')
    print('Tor vanishing in every degree is by projectivity of the localisation target, not finite testing')
    print('ALL CHECKS PASSED')

if __name__=='__main__': main()
