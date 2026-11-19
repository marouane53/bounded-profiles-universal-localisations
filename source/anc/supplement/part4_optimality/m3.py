"""Exact transfer and HH^3 coboundary test for finite thin dg algebras.

API: analyse(DGA) returns H^0, H^{-1}, transferred m3 and a coboundary test.
CLI: python3 m3.py algebra.json
JSON fields: names, degrees, unit, differential (row matrices), products.
products consists of records [i,j,output-vector]. Rational entries are strings
or integers. This tests formality only when cohomology lies in [-1,0].
It does not test homologicality of a proposed ring localisation.
"""
from itertools import product
import json,sys
from exact import *


def analyse(b):
    if any(d not in (-1,0,1) for d in b.degrees):
        raise ValueError('This implementation accepts a strict [-1,1] complex only')
    idx={d:[i for i,x in enumerate(b.degrees) if x==d] for d in (-1,0,1)}
    def emb(v,d):
        w=vec(b.n)
        for i,c in zip(idx[d],v):w[i]=c
        return w
    dm=[[b.d[i][j] for j in idx[-1]] for i in idx[0]]
    d0=[[b.d[i][j] for j in idx[0]] for i in idx[1]]
    h1=len(idx[1])-rank(d0,len(idx[0]))
    if h1:raise ValueError(f'Not a thin connective model: dim H^1={h1}')
    us=[emb(v,-1) for v in nullspace(dm,len(idx[-1]))]
    cm=[];bd=[]
    for j in idx[-1]:
        x=b.basis(j);y=matvec(b.d,x)
        if rank(hstack(bd+[y],b.n))>len(bd):cm.append(x);bd.append(y)
    zs=[emb(v,0) for v in nullspace(d0,len(idx[0]))]
    lifts=[]
    for x in [b.unit]+zs:
        if rank(hstack(bd+lifts+[x],b.n))>len(bd)+len(lifts):lifts.append(x)
    if not lifts:raise ValueError('The dg algebra is acyclic')
    c0=[]
    for j in idx[0]:
        x=b.basis(j)
        if rank(hstack(bd+lifts+c0+[x],b.n))>len(bd)+len(lifts)+len(c0):c0.append(x)
    basis0=bd+lifts+c0; basisminus=us+cm
    def ph0(x):return coords(basis0,x)[len(bd):len(bd)+len(lifts)]
    def pu(x):return coords(basisminus,x)[:len(us)]
    def h(x):return lincomb(cm,coords(basis0,x)[:len(bd)],b.n)
    g=len(lifts);u=len(us)
    mult=[[ph0(b.mul(x,y)) for y in lifts] for x in lifts]
    left=[[pu(b.mul(x,y)) for y in us] for x in lifts]
    right=[[pu(b.mul(y,x)) for y in us] for x in lifts]
    mu={}
    for a,c,d in product(range(g),repeat=3):
        mu[a,c,d]=pu(sub(b.mul(lifts[a],h(b.mul(lifts[c],lifts[d]))),
                          b.mul(h(b.mul(lifts[a],lifts[c])),lifts[d])))
    def l(a,v):return lincomb(left[a],v,u)
    def r(v,a):return lincomb(right[a],v,u)
    for a,c,d,e in product(range(g),repeat=4):
        z=l(a,mu[c,d,e])
        for j,t in enumerate(mult[a][c]):z=sub(z,scale(t,mu[j,d,e]))
        for j,t in enumerate(mult[c][d]):z=add(z,scale(t,mu[a,j,e]))
        for j,t in enumerate(mult[d][e]):z=sub(z,scale(t,mu[a,c,j]))
        z=add(z,r(mu[a,c,d],e))
        assert not any(z),('Hochschild cocycle failure',a,c,d,e,z)
    # Normalised cochains: lift 0 is the unit. The transferred cocycle is
    # normalised because the contraction annihilates the cohomology lifts.
    assert lifts[0]==b.unit
    for a,c,d in product(range(g),repeat=3):
        if 0 in (a,c,d):assert not any(mu[a,c,d])
    unknowns=list(product(range(1,g),range(1,g),range(u)))
    equations=list(product(range(1,g),range(1,g),range(1,g),range(u)))
    cols=[]
    for x,y,z in unknowns:
        def co(a,c):return vec(u,z) if (a,c)==(x,y) else vec(u)
        vals={}
        for a,c,d in product(range(1,g),repeat=3):
            v=l(a,co(c,d))
            for j,t in enumerate(mult[a][c]):v=sub(v,scale(t,co(j,d)))
            for j,t in enumerate(mult[c][d]):v=add(v,scale(t,co(a,j)))
            v=sub(v,r(co(a,c),d));vals[a,c,d]=v
        cols.append([vals[a,c,d][z] for a,c,d,z in equations])
    mat=hstack(cols,len(equations)); target=[mu[a,c,d][z] for a,c,d,z in equations]
    solution=solve(mat,target,len(unknowns))
    return dict(gamma_dim=g,u_dim=u,h1=h1,lifts=lifts,u_basis=us,
                product=mult,left=left,right=right,m3=mu,
                coboundary_rank=rank(mat,len(unknowns)),
                augmented_rank=rank([row+[v] for row,v in zip(mat,target)]),
                formal=solution is not None,coboundary_solution=solution)


def report(result):
    print('dim H^-1 =',result['u_dim'],'; dim H^0 =',result['gamma_dim'],'; dim H^1 =',result['h1'])
    nz={str(key):[str(x) for x in value] for key,value in result['m3'].items() if any(value)}
    print('Nonzero m3 entries in the chosen cohomology basis:',nz)
    print('Hochschild cocycle identity: PASS')
    print('rank(delta:C2->C3) =',result['coboundary_rank'])
    print('rank([delta | m3]) =',result['augmented_rank'])
    print('m3 is a Hochschild coboundary:',result['formal'])

if __name__=='__main__':
    if len(sys.argv)!=2:raise SystemExit('Usage: python3 m3.py algebra.json')
    with open(sys.argv[1]) as f:q=json.load(f)
    b=DGA(q['names'],q['degrees'],{(i,j):v for i,j,v in q['products']},q['unit'],q['differential'])
    b.check_dga();report(analyse(b))
