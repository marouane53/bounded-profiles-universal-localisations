#!/usr/bin/env python3
"""Exact finite certificates for the connective-realisation obstruction.

Arithmetic: fractions over Q, with separate prime-field rank checks.
This script does NOT purport to quantify over all rings/gluings. Those
nonexistence statements are proved in the accompanying mathematical text.
The helpers exact.py, m3.py and verify_nonformal.py are unmodified copies
of the user's supplied files. No finite range of Tor is used as a proof.
"""
from itertools import product
from fractions import Fraction as F
from pathlib import Path
import json
from exact import (Algebra, DGA, vec, add, sub, scale, lincomb, hstack,
                   matvec, matmul, rank, nullspace, independent)
from verify_nonformal import (example_dga, context_ring, endomorphism_dga)


def sign(n):
    return -1 if n % 2 else 1


def cohomology(degrees, differential):
    """Compute every cohomology degree of a finite complex, including zeros."""
    ans = {}
    for q in range(min(degrees), max(degrees) + 1):
        ids = [i for i, x in enumerate(degrees) if x == q]
        prev = [i for i, x in enumerate(degrees) if x == q - 1]
        nxt = [i for i, x in enumerate(degrees) if x == q + 1]
        incoming = [[differential[i][j] for j in prev] for i in ids]
        outgoing = [[differential[i][j] for j in ids] for i in nxt]
        ans[q] = len(ids) - rank(incoming, len(prev)) - rank(outgoing, len(ids))
    return ans


def nonzero_h(b):
    return {q: h for q, h in cohomology(b.degrees, b.d).items() if h}


def cohomology_lifts(b, q):
    ids = [i for i, d in enumerate(b.degrees) if d == q]
    nxt = [i for i, d in enumerate(b.degrees) if d == q + 1]
    local = [[b.d[i][j] for j in ids] for i in nxt]
    cycles = []
    for x in nullspace(local, len(ids)):
        out = vec(b.n)
        for i, c in zip(ids, x): out[i] = c
        cycles.append(out)
    boundaries = independent([matvec(b.d, b.basis(i))
                              for i, d in enumerate(b.degrees) if d == q - 1], b.n)
    basis = list(boundaries)
    lifts = []
    for x in cycles:
        if rank(hstack(basis + [x], b.n)) > len(basis):
            basis.append(x); lifts.append(x)
    return lifts, boundaries


def square_zero_generator(degree, name):
    return DGA(['1', name], [0, degree],
               {(0, 0): {0: 1}, (0, 1): {1: 1}, (1, 0): {1: 1}},
               [1, 0], [[0, 0], [0, 0]])


def exterior_two():
    tab = {}
    for i in range(4): tab[0, i] = {i: 1}; tab[i, 0] = {i: 1}
    tab[1, 2] = {3: 1}; tab[2, 1] = {3: -1}
    return DGA(['1', 'u', 'v', 'uv'], [0, -1, -1, -2], tab,
               [1, 0, 0, 0], [[0] * 4 for _ in range(4)])


def tensor_dga(a, b):
    tags = list(product(range(a.n), range(b.n)))
    ix = {p: i for i, p in enumerate(tags)}
    n = len(tags)
    tab = {}
    for f, (i, j) in enumerate(tags):
        for g, (r, s) in enumerate(tags):
            out = vec(n)
            for l, c in enumerate(a.mul(a.basis(i), a.basis(r))):
                if not c: continue
                for m, z in enumerate(b.mul(b.basis(j), b.basis(s))):
                    if z: out[ix[l, m]] += sign(b.degrees[j] * a.degrees[r]) * c * z
            if any(out): tab[f, g] = out
    dcols = []
    for i, j in tags:
        out = vec(n)
        for l, c in enumerate(matvec(a.d, a.basis(i))):
            if c: out[ix[l, j]] += c
        for m, c in enumerate(matvec(b.d, b.basis(j))):
            if c: out[ix[i, m]] += sign(a.degrees[i]) * c
        dcols.append(out)
    unit = [a.unit[i] * b.unit[j] for i, j in tags]
    d = DGA([a.names[i] + '*' + b.names[j] for i, j in tags],
            [a.degrees[i] + b.degrees[j] for i, j in tags], tab, unit, hstack(dcols, n))
    d.tags = tags; d.lookup = ix
    return d


def fibre_of_degree_zero_inclusion(d):
    """Cone(D^0 -> D)[-1]: D^{q-1} plus (D^0)^q, d(z,a)=(-dz+a,0)."""
    ai = [i for i, q in enumerate(d.degrees) if q == 0]
    degrees = [q + 1 for q in d.degrees] + [0] * len(ai)
    n = len(degrees)
    dm = [[F(0)] * n for _ in range(n)]
    for i in range(d.n):
        for j in range(d.n): dm[i][j] = -d.d[i][j]
    for j, i in enumerate(ai): dm[i][d.n + j] = 1
    assert not any(x for row in matmul(dm, dm) for x in row)
    return {q: n for q, n in cohomology(degrees, dm).items() if n}


def check_dg_map(a, b, cols):
    assert len(cols) == a.n
    phi = lambda x: lincomb(cols, x, b.n)
    assert phi(a.unit) == b.unit
    for i in range(a.n):
        assert all(not c or b.degrees[j] == a.degrees[i] for j, c in enumerate(cols[i]))
        assert phi(matvec(a.d, a.basis(i))) == matvec(b.d, cols[i])
        for j in range(a.n):
            assert phi(a.mul(a.basis(i), a.basis(j))) == b.mul(cols[i], cols[j])
    assert nonzero_h(a) == nonzero_h(b), (nonzero_h(a), nonzero_h(b))
    for q, dim in nonzero_h(a).items():
        lifts, _ = cohomology_lifts(a, q)
        _, boundaries = cohomology_lifts(b, q)
        assert rank(hstack(boundaries + [phi(x) for x in lifts], b.n)) == len(boundaries) + dim


def radical_square_zero_cycle(n):
    names = [f'e{i+1}' for i in range(n)] + [f'a{i+1}' for i in range(n)]
    tab = {}
    for i in range(n):
        tab[i, i] = {i: 1}
        tab[i, n+i] = {n+i: 1}
        tab[n+i, (i+1) % n] = {n+i: 1}
    r = Algebra(names, tab, [1]*n + [0]*n)
    r.starts = list(range(n)) + list(range(n))
    r.ends = list(range(n)) + [(i+1) % n for i in range(n)]
    r.idempotents = [r.basis(i) for i in range(n)]
    return r


def ordinary_tensor_ring(a):
    tags = list(product(range(a.n), repeat=2)); ix = {p: i for i, p in enumerate(tags)}
    n = len(tags); tab = {}
    for f, (i, j) in enumerate(tags):
        for g, (r, s) in enumerate(tags):
            out = vec(n)
            for l, c in enumerate(a.mul(a.basis(i), a.basis(r))):
                if c:
                    for m, z in enumerate(a.mul(a.basis(j), a.basis(s))):
                        if z: out[ix[l, m]] += c*z
            if any(out): tab[f, g] = out
    r = Algebra([a.names[i]+'|'+a.names[j] for i,j in tags], tab,
                [a.unit[i]*a.unit[j] for i,j in tags])
    r.starts = [2*a.starts[i]+a.starts[j] for i,j in tags]
    r.ends = [2*a.ends[i]+a.ends[j] for i,j in tags]
    r.idempotents = []
    for p,q in product(range(2), repeat=2):
        r.idempotents.append([a.idempotents[p][i]*a.idempotents[q][j] for i,j in tags])
    r.tags=tags; r.lookup=ix
    return r


def underlying_complex(r, summands, differential):
    tags = [(j, i) for j,(v,q) in enumerate(summands) for i in range(r.n) if r.ends[i] == v]
    ix = {x:i for i,x in enumerate(tags)}; n=len(tags)
    deg = [summands[j][1] for j,i in tags]
    dm = [[F(0)]*n for _ in range(n)]
    for c,(j,i) in enumerate(tags):
        for (src,tgt),x in differential.items():
            if src != j: continue
            for h,z in enumerate(r.mul(r.basis(i),x)):
                if z: dm[ix[tgt,h]][c] += z
    assert not any(x for row in matmul(dm,dm) for x in row)
    return tags,deg,dm


def check_components(r, summands, differential, survivor, expected):
    tags,deg,dm=underlying_complex(r,summands,differential)
    ans={}
    for v in range(len(r.idempotents)):
        ids=[j for j,(_,i) in enumerate(tags) if r.starts[i]==v]
        ds=[deg[j] for j in ids]
        dd=[[dm[i][j] for j in ids] for i in ids]
        h={q:m for q,m in cohomology(ds,dd).items() if m}
        ans[v]=h
        assert h == (expected if v==survivor else {})
    return ans


def hom_to_projective(r, summands, differential, target_vertex):
    # Use the actual mapping complex, not End dimensions inferred from theory.
    both=summands+[(target_vertex,0)]
    e=endomorphism_dga(r,both,differential)
    ids=[i for i,(src,tgt,q) in enumerate(e.map_tags) if src<len(summands) and tgt==len(summands)]
    degrees=[e.degrees[i] for i in ids]
    dm=[[e.d[i][j] for j in ids] for i in ids]
    return {q:n for q,n in cohomology(degrees,dm).items() if n}


def dump_ring(r, name):
    path = Path(__file__).parent/(name+'_multiplication.json')
    data={'basis':r.names, 'unit':[str(x) for x in r.unit],
          'products':[[i,j,[str(x) for x in v]] for (i,j),v in sorted(r.table.items()) if any(v)],
          'omitted_products':'zero', 'starts':r.starts, 'ends':r.ends}
    path.write_text(json.dumps(data,indent=2)+'\n')


def rank_mod_p(matrix, p):
    a=[[int(x.numerator)*pow(int(x.denominator),-1,p)%p for x in row] for row in matrix]
    if not a: return 0
    row=0
    for c in range(len(a[0])):
        pivot=next((i for i in range(row,len(a)) if a[i][c]),None)
        if pivot is None: continue
        a[row],a[pivot]=a[pivot],a[row]
        inv=pow(a[row][c],-1,p); a[row]=[(x*inv)%p for x in a[row]]
        for i in range(len(a)):
            if i!=row:
                coeff=a[i][c]; a[i]=[(x-coeff*y)%p for x,y in zip(a[i],a[row])]
        row+=1
        if row==len(a):break
    return row


def cohom_mod_p(b,p):
    out={}
    for q in range(min(b.degrees),max(b.degrees)+1):
        ids=[i for i,d in enumerate(b.degrees) if d==q]
        prev=[i for i,d in enumerate(b.degrees) if d==q-1]
        nxt=[i for i,d in enumerate(b.degrees) if d==q+1]
        hi=len(ids)-rank_mod_p([[b.d[i][j] for j in prev] for i in ids],p)-rank_mod_p([[b.d[i][j] for j in ids] for i in nxt],p)
        if hi:out[q]=hi
    return out


def main():
    thin=square_zero_generator(-1,'eta'); sphere=square_zero_generator(-2,'u'); ext=exterior_two()
    d0=example_dga(); nf=tensor_dga(d0,thin)
    examples=[('thin control',thin),('sphere degree -2',sphere),('exterior two',ext),('nonformal tensor',nf)]
    print('Exact connective-realisation certificates over Q')
    for label,d in examples:
        d.check_dga()
        print(label+': H = '+str(nonzero_h(d))+'; fibre H = '+str(fibre_of_degree_zero_inclusion(d)))
    assert nonzero_h(nf)=={-2:1,-1:3,0:2}
    assert fibre_of_degree_zero_inclusion(thin)=={0:1}
    assert fibre_of_degree_zero_inclusion(sphere)=={-1:1}
    assert fibre_of_degree_zero_inclusion(ext)=={-1:1,0:2}
    assert fibre_of_degree_zero_inclusion(nf)=={-1:1,0:6}
    # Nonformal Massey witness. The indeterminacy is computed, not stipulated.
    t=nf.basis(nf.lookup[1,0]); r=nf.basis(nf.lookup[7,0]); w=nf.basis(nf.lookup[8,0])
    assert matvec(nf.d,r)==nf.mul(t,t)
    witness=sub(nf.mul(t,r),nf.mul(r,t)); assert witness==w
    lifts,boundaries=cohomology_lifts(nf,-1)
    indet=[nf.mul(t,x) for x in lifts]+[nf.mul(x,t) for x in lifts]
    base_rank=rank(hstack(boundaries,nf.n))
    ind_rank=rank(hstack(boundaries+indet,nf.n))
    wit_rank=rank(hstack(boundaries+indet+[w],nf.n))
    assert ind_rank-base_rank==1 and wit_rank==ind_rank+1
    print('Nonformal tensor: <t,t,t> represented by w; indeterminacy dimension 1; 0 excluded: PASS')
    # Frobenius/graded-symmetry certificates for the Morita-invariant obstruction.
    for label,d,top in [('sphere',sphere,1),('exterior',ext,3)]:
        gram=[[d.mul(d.basis(i),d.basis(j))[top] for j in range(d.n)] for i in range(d.n)]
        assert rank(gram)==d.n
        for i,j in product(range(d.n),repeat=2):
            assert gram[i][j]==sign(d.degrees[i]*d.degrees[j])*gram[j][i]
        print(label+': degree -2 graded-symmetric Frobenius pairing rank '+str(rank(gram))+': PASS')
    # Correct the two-term template on both formal and nonformal controls.
    for label,d in [('formal',thin),('nonformal',d0)]:
        rr=context_ring(d); sums=[(0,-1),(1,0)]; diff={(0,1):rr.s}
        hp=hom_to_projective(rr,sums,diff,0); he=hom_to_projective(rr,sums,diff,1)
        assert hp=={q+1:n for q,n in nonzero_h(d).items()} and not he
        print(label+' two-term: H Theta(Rp) = '+str(hp)+' = H(D[-1]); H Theta(Re) = {}: PASS')
    # Explicit Schur realisation of the degree -2 sphere.
    a=radical_square_zero_cycle(3); a.check()
    sums=[(1,-2),(2,-1),(0,0)]; diff={(0,1):a.basis(4),(1,2):a.basis(5)}
    end=endomorphism_dga(a,sums,diff); end.check_dga()
    u=vec(end.n);u[end.map_lookup[2,0,3]]=1
    check_dg_map(sphere,end,[end.unit,u])
    hc=cohomology(end.degrees,end.d)
    comp=check_components(a,sums,diff,0,nonzero_h(sphere))
    print('Sphere Schur ring: dim R = 6; dim eRe = 3; End(T)^op H by every degree = '+str(hc))
    print('Sphere strict dg map and idempotent-component cohomology = '+str(comp)+': PASS')
    dump_ring(a,'sphere_schur')
    # Tensor of the two-term formal realisation, for the exterior algebra.
    c=context_ring(thin); rr=ordinary_tensor_ring(c);rr.check()
    sums=[(0,-2),(1,-1),(2,-1),(3,0)]
    ix=rr.lookup
    p=next(i for i,x in enumerate(c.tags) if x==('a',0))
    e=next(i for i,x in enumerate(c.tags) if x==('b',0))
    s=next(i for i,x in enumerate(c.tags) if x==('v',0))
    m=next(i for i,x in enumerate(c.tags) if x==('m',1))
    diff={(0,2):rr.basis(ix[s,p]), (0,1):scale(-1,rr.basis(ix[p,s])),
          (1,3):rr.basis(ix[s,e]), (2,3):rr.basis(ix[e,s])}
    ee=endomorphism_dga(rr,sums,diff);ee.check_dga()
    U=vec(ee.n); V=vec(ee.n)
    U[ee.map_lookup[2,0,ix[m,p]]]=1;U[ee.map_lookup[3,1,ix[m,e]]]=1
    V[ee.map_lookup[1,0,ix[p,m]]]=-1;V[ee.map_lookup[3,2,ix[e,m]]]=1
    UV=ee.mul(U,V)
    check_dg_map(ext,ee,[ee.unit,U,V,UV])
    comp=check_components(rr,sums,diff,3,nonzero_h(ext))
    corner=sum(1 for x,y in zip(rr.starts,rr.ends) if x!=3 and y!=3)
    print('Exterior Schur ring: dim R = 16; dim hRh = '+str(corner)+'; End(T)^op H by every degree = '+str(cohomology(ee.degrees,ee.d)))
    print('Exterior strict dg map, uv = -vu, and component cohomology = '+str(comp)+': PASS')
    dump_ring(rr,'exterior_schur')
    for p in [2,3,101]:
        assert cohom_mod_p(end,p)==nonzero_h(sphere)
        assert cohom_mod_p(ee,p)==nonzero_h(ext)
        assert cohom_mod_p(nf,p)==nonzero_h(nf)
    print('Independent rank reductions over F_2, F_3, F_101: PASS')
    print('All finite complexes checked in their ENTIRE degree ranges; no truncated Tor test used.')
    print('ALL CHECKS PASSED')

if __name__=='__main__': main()
