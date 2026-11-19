"""Cross-check noncommutative rewriting, open-chart coordinates and finite resolutions."""
from collections import Counter
from functools import cache
from itertools import product
from check_lib import rank, sparse_rank, tidy, scaled_add, qpower, ident, mul, require_asserts


def normal_form_engine(n,Q,p):
    @cache
    def normalize(w,last=False):
        sites=[i for i in range(len(w)-n) if w[i:i+n]==(0,)*n and w[i+n]]
        if not sites:return {w:1}
        i=sites[-1] if last else sites[0]
        out=Counter()
        for a,c in enumerate(Q[w[i+n]-1],1):
            if c:
                v=w[:i]+(a,)+(0,)*n+w[i+n+1:]
                scaled_add(out,normalize(v,last),c)
        return tidy(out,p)
    return normalize


def normal_checker(n,Q,p,depth,lifts=False):
    r=len(Q)
    normalize=normal_form_engine(n,Q,p)
    powers=[qpower(Q,i,p) for i in range(depth+2)]
    def chart(w):
        # Degree-zero word times x^len(w) in the skew Laurent algebra.
        out={():1}
        for j,a in enumerate(w):
            if a:
                nxt=Counter()
                for b,c in enumerate(powers[j//n][a-1]):
                    if c:
                        for v,alpha in out.items():
                            nxt[v+((b,j%n),)]+=alpha*c
                out=tidy(nxt,p)
        return out
    bases=[]
    for d in range(depth+1):
        basis=[]
        for w in product(range(r+1),repeat=d):
            first=normalize(w)
            assert first==normalize(w,True)
            if first=={w:1}:basis.append(w)
            actual=Counter()
            for v,k in first.items():scaled_add(actual,chart(v),k)
            assert tidy(actual,p)==chart(w)
        bases.append(basis)
        recursive=1+sum(r*len(bases[d-a]) for a in range(1,min(n,d)+1))
        assert len(basis)==recursive
    if not lifts:return
    # Two-dimensional genuinely matrix-valued lifts; no commutation imposed.
    ops=[tuple(tuple((1+(i+1)*(a+1)+(j+2)*(b+1)+int(a==b)*(j+i)) % p
                    for b in range(2)) for a in range(2))
         for i,j in product(range(r),range(n))]
    def action(w):
        mat=ident(2)
        for j,a in enumerate(w):
            if not a:continue
            mix=powers[j//n][a-1]
            nxt=tuple(tuple(sum(co*ops[b*n+j%n][s][t] for b,co in enumerate(mix)) % p
                            for t in range(2)) for s in range(2))
            mat=mul(mat,nxt,p)
        return mat
    for d,basis in enumerate(bases):
        columns=[]
        for a in range(min(d,n)):
            for i in range(r):
                for suffix in bases[d-a-1]:
                    for v in range(2):
                        col=Counter({(v,(0,)*a+(i+1,)+suffix):1})
                        for j,c in enumerate(ops[i*n+a][v]):
                            for w,alpha in normalize((0,)*(a+1)+suffix).items():
                                col[j,w]-=c*alpha
                        col=tidy(col,p)
                        image=[sum(c*action(w)[v][j] for (v,w),c in col.items()) % p for j in range(2)]
                        assert image==[0,0]
                        columns.append(col)
        assert len(columns)==sparse_rank(columns,p)==2*len(basis)-2


def main():
    runs=0
    for p in (2,3):
        for coeffs in product(range(p),repeat=4):
            Q=(coeffs[:2],coeffs[2:])
            if rank(Q,p)!=2:continue
            for n in (1,2,3,4):
                normal_checker(n,Q,p,4)
                runs+=1
    for p in (2,3,5):
        for n in (1,2,3,4):
            normal_checker(n,((1,1),(0,1)),p,6,True)
            normal_checker(n,((0,1),(1,0)),p,5,True)
        normal_checker(3,((1,1,0),(0,1,1),(0,0,1)),p,5,True)
    print(f'{runs} exhaustive invertible matrix/period cases over F2,F3, all words through length 4.')
    print('27 additional rank-two module resolutions and rewrite checks over F2,F3,F5, periods 1–4 and up to 3 z-generators.')
    print('Included unipotent twists of characteristic-divisible order and odd normal elements in odd characteristic.')
    # Inseparable simple over k=F2(t): alpha²=t has no root in k by the t-adic valuation.
    # Compute the 2x2 companion determinant symbolically using bit polynomials over F2.
    # det(alpha I - C) = alpha²+t = 0; det((alpha+1)I-C)=1 in the same field.
    # Packed polynomials represent F2[t]; pairs represent a+b*alpha.
    def times(a,b):
        out=0
        while b:
            if b&1:out^=a
            a<<=1
            b>>=1
        return out
    def extension_square(a,b):
        return times(a,a) ^ (times(b,b)<<1)
    assert extension_square(0,1)^2==0 and extension_square(1,1)^2==1
    print('Companion separator checked in F2(t)[alpha]/(alpha²−t); the distinct translate gives determinant 1.')


if __name__=='__main__':
    require_asserts()
    main()
