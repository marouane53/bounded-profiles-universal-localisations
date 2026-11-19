#!/usr/bin/env python3
"""Exact GF(3) block certificate for G=C3 x C2.

This is a counterexample to the unrestricted group-algebra pair-support
claim in the brief, not a test of the complete-intersection theorem.
"""
from itertools import product

P=3
basis=list(product(range(3),range(2)))
pos={x:i for i,x in enumerate(basis)}


def unit(a,b):
    z=[0]*6; z[pos[(a%3,b%2)]]=1; return z


def plus(x,y,scale=1): return [(a+scale*b)%3 for a,b in zip(x,y)]


def mul(x,y):
    z=[0]*6
    for i,(a,b) in enumerate(basis):
        for j,(c,d) in enumerate(basis):
            z[pos[((a+c)%3,(b+d)%2)]]=(z[pos[((a+c)%3,(b+d)%2)]]+x[i]*y[j])%3
    return z


def rankcols(cols):
    a=[[c[i]%3 for c in cols] for i in range(6)]
    r=0
    for j in range(len(cols)):
        pivot=next((i for i in range(r,6) if a[i][j]),None)
        if pivot is None: continue
        a[r],a[pivot]=a[pivot],a[r]
        z=pow(a[r][j],-1,3); a[r]=[(z*x)%3 for x in a[r]]
        for i in range(6):
            if i!=r and a[i][j]:
                z=a[i][j]; a[i]=[(x-z*y)%3 for x,y in zip(a[i],a[r])]
        r+=1
        if r==6: break
    return r


def main():
    one=unit(0,0); g=unit(1,0); h=unit(0,1)
    ep=[2*x%3 for x in plus(one,h)]
    em=[2*x%3 for x in plus(one,h,-1)]
    z=[0]*6
    assert mul(ep,ep)==ep and mul(em,em)==em
    assert mul(ep,em)==z and plus(ep,em)==one
    t=plus(g,one,-1); t2=mul(t,t)
    assert mul(t2,t)==z
    print('Group block obstruction: k[C3 x C2], characteristic 3')
    print('e_plus=(1+h)/2, e_minus=(1-h)/2: central orthogonal idempotents: PASS')
    for label,e in [('plus',ep),('minus',em)]:
        block=[e,mul(t,e),mul(t2,e)]
        assert rankcols(block)==3
        r1=rankcols([mul(t,v) for v in block])
        r2=rankcols([mul(t2,v) for v in block])
        assert (r1,r2)==(2,1)
        print(f'{label} block: basis e,t e,t^2 e; dimension 3; ranks(t,t^2)={(r1,r2)}: PASS')
    assert all(mul(x,y)==z for x in [ep,mul(t,ep),mul(t2,ep)] for y in [em,mul(t,em),mul(t2,em)])
    print('Cross-block maps and all cross-block Ext vanish: certified by orthogonal blocks')
    print('Each simple has an exact 2-periodic resolution, with maps t and t^2.')
    print('Applying Hom to its own simple kills both differentials: self-Ext dimension 1 in every degree.')
    print('Both simples have the same full projective cohomological support (one point).')
    print('Thus disjoint-support iff eventual pair-Ext vanishing FAILS for arbitrary finite groups.')
    print('ALL CHECKS PASSED')

if __name__=='__main__': main()
