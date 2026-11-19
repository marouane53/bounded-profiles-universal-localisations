#!/usr/bin/env python3
"""Nonformal thin kernel: exact ring, bimodule, cone, and HH^3 certificates.
The all-degree stable-flatness proof is the two projective decompositions,
not a finite Tor truncation.  All computations are over Q.
"""
from itertools import product
from exact import *
from m3 import analyse,report


def example_dga():
    # A basis 1,t,a,b,c; M basis z,l,r,w.  H^-1 is k*w.
    names=['1','t','a','b','c','z','l','r','w'];n=len(names)
    tab={}
    for j in range(n):tab[0,j]={j:1};tab[j,0]={j:1}
    tab.update({(1,1):{4:1},(1,2):{3:1},(2,1):{4:1},
                (1,5):{6:1},(1,7):{8:1},(5,1):{7:1},(6,1):{8:1}})
    d=[vec(n) for _ in range(n)]
    for src,tgt in [(5,2),(6,3),(7,4)]:d[tgt][src]=1
    return DGA(names,[0]*5+[-1]*4,tab,vec(n,0),d)


def context_ring(dga):
    ai=[i for i,n in enumerate(dga.degrees) if n==0]
    mi=[i for i,n in enumerate(dga.degrees) if n==-1]
    assert set(dga.degrees)<=set([-1,0])
    tags=[(corner,i) for corner in ('a','v','m','b') for i in (mi if corner=='m' else ai)]
    lookup={x:i for i,x in enumerate(tags)};n=len(tags)
    def put(corner,value):
        out=vec(n)
        for i,c in enumerate(value):
            if c:out[lookup[corner,i]]=c
        return out
    tab={}
    for r,(x,i) in enumerate(tags):
        for s,(y,j) in enumerate(tags):
            corner=None;v=vec(dga.n)
            if (x,y) in [('a','a'),('a','v'),('v','b'),('m','a'),('b','m'),('b','b')]:
                corner={('a','a'):'a',('a','v'):'v',('v','b'):'v',('m','a'):'m',('b','m'):'m',('b','b'):'b'}[x,y]
                v=dga.mul(dga.basis(i),dga.basis(j))
            elif (x,y)==('v','m'):
                corner='a';v=dga.mul(dga.basis(i),matvec(dga.d,dga.basis(j)))
            elif (x,y)==('m','v'):
                corner='b';v=dga.mul(matvec(dga.d,dga.basis(i)),dga.basis(j))
            if corner and any(v):tab[r,s]=put(corner,v)
    unit=add(put('a',dga.unit),put('b',dga.unit))
    r=Algebra([x+':'+dga.names[i] for x,i in tags],tab,unit)
    r.starts=[0 if x in ('a','v') else 1 for x,i in tags]
    r.ends=[0 if x in ('a','m') else 1 for x,i in tags]
    r.idempotents=[put('a',dga.unit),put('b',dga.unit)]
    r.s=put('v',dga.unit)
    r.tags=tags;r.lookup=lookup
    # Matrix algebra M_2(A), with A = D^0.
    stags=[(h,l,i) for h in (0,1) for l in (0,1) for i in ai]
    si={v:j for j,v in enumerate(stags)};sn=len(stags)
    stab={}
    for q,(h,l,i) in enumerate(stags):
        for t,(h2,l2,j) in enumerate(stags):
            if l!=h2:continue
            v=dga.mul(dga.basis(i),dga.basis(j));out=vec(sn)
            for k,c in enumerate(v):
                if c:out[si[h,l2,k]]=c
            if any(out):stab[q,t]=out
    su=vec(sn)
    for h in (0,1):
        for i,c in enumerate(dga.unit):
            if c:su[si[h,h,i]]=c
    s=Algebra([f'{h}{l}:{dga.names[i]}' for h,l,i in stags],stab,su)
    lam=[]
    for corner,i in tags:
        out=vec(sn)
        if corner=='m':
            for j,c in enumerate(matvec(dga.d,dga.basis(i))):
                if c:out[si[1,0,j]]=c
        else:out[si[{'a':0,'v':0,'b':1}[corner],{'a':0,'v':1,'b':1}[corner],i]]=1
        lam.append(out)
    r.lambda_columns=lam;r.target=s;r.target_tags=stags
    return r


def endomorphism_dga(r, summands, differential):
    """End(T)^op for a finite complex of projectives R e_v.

    r must have a Peirce-adapted basis, with starts, ends and idempotents.
    summands is [(vertex,cohomological degree),...].
    differential maps (source summand,target summand) to a coefficient in R.
    This supports matrices and direct sums, not just one indecomposable map.
    """
    tags=[]
    for src,(v,ds) in enumerate(summands):
        for tgt,(w,dt) in enumerate(summands):
            for q in range(r.n):
                if r.starts[q]==v and r.ends[q]==w:tags.append((src,tgt,q))
    ix={x:j for j,x in enumerate(tags)};n=len(tags)
    deg=[summands[t][1]-summands[s][1] for s,t,q in tags]
    def compose_after(f,g):
        # g after f, underlying ungraded composition, for right coefficients.
        a,b,q=tags[f];c,d,t=tags[g];out=vec(n)
        if b!=c:return out
        for k,z in enumerate(r.mul(r.basis(q),r.basis(t))):
            if z:out[ix[a,d,k]]=z
        return out
    tab={}
    for f,g in product(range(n),repeat=2):
        v=scale((-1 if (deg[f]*deg[g]) % 2 else 1),compose_after(f,g))
        if any(v):tab[f,g]=v
    u=vec(n)
    for j,(v,d) in enumerate(summands):
        for q,c in enumerate(r.idempotents[v]):
            if c:u[ix[j,j,q]]=c
    dt=vec(n)
    for (a,b),v in differential.items():
        for q,c in enumerate(v):
            if c:dt[ix[a,b,q]]=c
    dc=[]
    for f in range(n):
        out=vec(n)
        for g,c in enumerate(dt):
            if c:out=add(out,scale(c,sub(compose_after(f,g),scale((-1 if deg[f] % 2 else 1),compose_after(g,f)))))
        dc.append(out)
    out=DGA([f'{a}>{b}:{r.names[q]}' for a,b,q in tags],deg,tab,u,hstack(dc,n))
    out.map_tags=tags;out.map_lookup=ix
    return out


def verify(d):
    d.check_dga()
    r=context_ring(d);s=r.target
    r.check();s.check()
    la=lambda x:lincomb(r.lambda_columns,x,s.n)
    for i,j in product(range(r.n),repeat=2):
        assert la(r.mul(r.basis(i),r.basis(j)))==s.mul(la(r.basis(i)),la(r.basis(j)))
    assert la(r.unit)==s.unit
    # Row and column identifications with pR and Re, on every basis action.
    rowids=[i for i in range(r.n) if r.starts[i]==0]
    colids=[i for i in range(r.n) if r.ends[i]==1]
    st=r.target_tags; sti={v:i for i,v in enumerate(st)}
    def row_map(x,row):
        out=vec(s.n)
        for i,c in enumerate(x):
            if c:
                corner,q=r.tags[i];assert corner in ('a','v')
                out[sti[row,0 if corner=='a' else 1,q]]=c
        return out
    def col_map(x,col):
        out=vec(s.n)
        for i,c in enumerate(x):
            if c:
                corner,q=r.tags[i];assert corner in ('v','b')
                out[sti[0 if corner=='v' else 1,col,q]]=c
        return out
    for a in range(r.n):
        for x in rowids:
            for row in (0,1):assert row_map(r.mul(r.basis(x),r.basis(a)),row)==s.mul(row_map(r.basis(x),row),la(r.basis(a)))
        for x in colids:
            for col in (0,1):assert col_map(r.mul(r.basis(a),r.basis(x)),col)==s.mul(la(r.basis(a)),col_map(r.basis(x),col))
    b=endomorphism_dga(r,[(0,-1),(1,0)],{(0,1):r.s});b.check_dga()
    incl=[]
    for j,degree in enumerate(d.degrees):
        x=vec(b.n)
        if degree==0:
            x[b.map_lookup[0,0,r.lookup['a',j]]]=1
            x[b.map_lookup[1,1,r.lookup['b',j]]]=1
        else:x[b.map_lookup[1,0,r.lookup['m',j]]]=1
        incl.append(x)
    im=lambda x:lincomb(incl,x,b.n)
    for i in range(d.n):
        assert im(matvec(d.d,d.basis(i)))==matvec(b.d,im(d.basis(i)))
        for j in range(d.n):assert im(d.mul(d.basis(i),d.basis(j)))==b.mul(im(d.basis(i)),im(d.basis(j)))
    hb=analyse(b);hd=analyse(d)
    assert (hb['gamma_dim'],hb['u_dim'],hb['formal'])==(hd['gamma_dim'],hd['u_dim'],hd['formal'])
    return r,s,b,hb


def main():
    d=example_dga();r,s,b,h=verify(d)
    print('Nonformal crossed-module localisation over Q')
    print('dim A = 5; dim M = 4; dim R =',r.n,'; dim S =',s.n)
    print('Associativity, bimodule laws, Peiffer identity and dga Leibniz: PASS')
    print('lambda multiplicative on all basis pairs: PASS')
    print('S_R = (pR)^2 and _R S = (Re)^2, all basis actions: PASS')
    print('Strict dg map D -> End_R(T)^op, all products and differentials: PASS')
    print('Transferred structure on D:')
    report(analyse(d))
    print('Transferred structure computed independently from the presentation cone:')
    report(h)
    assert not h['formal']
    # A zero-differential comparison with the same H^0 and H^-1.
    tab={(0,0):{0:1},(0,1):{1:1},(1,0):{1:1},(0,2):{2:1},(2,0):{2:1}}
    formal=DGA(['1','t','u'],[0,0,-1],tab,[1,0,0],[[0]*3 for _ in range(3)])
    _,_,_,hf=verify(formal);assert hf['formal']
    print('Formal same-cohomology comparison: m3 = 0: PASS')
    print('ALL CHECKS PASSED')

if __name__=='__main__':main()
