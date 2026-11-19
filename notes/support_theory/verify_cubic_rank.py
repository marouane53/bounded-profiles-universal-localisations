"""Compute minimal and strict bar ranks of the cubic matrix family over prime fields."""
from collections import Counter
from functools import cache
from itertools import product
from check_lib import rank, sparse_rank, sign, tidy, scaled_add, weighted_words, require_asserts


def predicted(r,m,d,n):
    dims=[1]
    for i in range(1,n+1):
        value=(r+1)*dims[i-1]
        if i>=2: value+=(m-d)*dims[i-2]
        if i>=3: value-=d*dims[i-3]
        dims.append(value)
    return dims


def reduced_ranks(C,p,n):
    m,r=len(C),len(C[0])
    weights=(1,)*(r+1)+(2,)*m
    @cache
    def delta(w):
        out=Counter()
        parity=1
        for i,a in enumerate(w):
            if a<=r:
                parity=-parity
                continue
            for b,c in enumerate(C[a-r-1],1):
                if c:
                    out[w[:i]+(0,0,b)+w[i+1:]]+=parity*c
                    out[w[:i]+(b,0,0)+w[i+1:]]-=parity*c
        return tidy(out,p)
    values=[]
    before=0
    for i in range(n+1):
        words=weighted_words(weights,i)
        columns=[delta(w) for w in words]
        after=sparse_rank(columns,p)
        values.append(len(words)-before-after)
        before=after
    return values


def strict_ranks(C,p,n):
    # Original radical basis: a, b_i, h (degree zero), e, w_j (degree -1).
    # Each suspension contributes 1 to bar degree, plus the absolute internal degree.
    m,r=len(C),len(C[0])
    h,e=r+1,r+2
    weights=(1,)*(r+2)+(2,)*(m+1)
    def mult(a,b):
        if a==b==0:return {h:1}
        if a==e and 1<=b<=r: a,b=b,a
        if b==e and 1<=a<=r:
            return {e+j+1:row[a-1] for j,row in enumerate(C) if row[a-1]}
        return {}
    @cache
    def boundary(w):
        out=Counter()
        parity=1
        for i,a in enumerate(w):
            if a==e:
                out[w[:i]+(h,)+w[i+1:]]-=parity
            if i+1<len(w):
                for k,c in mult(a,w[i+1]).items():
                    out[w[:i]+(k,)+w[i+2:]]+=parity*sign(weights[a]-1)*c
            parity*=sign(weights[a])
        return tidy(out,p)
    ranks=[]
    for i in range(n+2):
        cols=[]
        for w in weighted_words(weights,i):
            col=boundary(w)
            sq=Counter()
            for b,c in col.items(): scaled_add(sq,boundary(b),c)
            assert not tidy(sq,p)
            cols.append(col)
        ranks.append(sparse_rank(cols,p))
    return [len(weighted_words(weights,i))-ranks[i]-ranks[i+1] for i in range(n+1)]


def main():
    grand=0
    for p,m,r,lim in [(2,1,3,5),(2,3,1,5),(2,2,2,6),(3,2,2,5),(5,2,2,4),
                       (2,2,3,4),(2,3,2,4),(2,3,3,4)]:
        counts=Counter()
        for entries in product(range(p),repeat=m*r):
            C=tuple(entries[i*r:(i+1)*r] for i in range(m))
            d=rank(C,p)
            observed=reduced_ranks(C,p,lim)
            assert observed==predicted(r,m,d,lim),(p,C,observed)
            counts[d]+=1
        grand+=sum(counts.values())
        print(f'F{p} {m}x{r}, all {sum(counts.values())} matrices through degree {lim}, rank histogram {dict(sorted(counts.items()))}.')
    strict_cases=[
        ((0,0),(0,0)),((1,1),(1,1)),((0,1),(1,0)),((1,2),(2,1)),
        ((1,1,1),(0,1,2)),((1,),(2,),(1,)),((1,0,1),(0,1,1),(1,1,0))]
    for p in (2,3,5):
        for C in strict_cases:
            out=strict_ranks(C,p,4)
            assert out==predicted(len(C[0]),len(C),rank(C,p),4)
    print(f'{grand} exhaustive matrix/field instances; also 21 original strict-bar calculations through degree 4, with d² checks.')


if __name__=='__main__':
    require_asserts()
    main()
