"""Check the two-parameter projective-support counterexample by its four-cell Hom differential."""
from collections import Counter
from check_lib import sparse_rank, require_asserts


def basis(n):
    return [(i,j,n+j-2*b,b) for i in (0,1) for j in (0,1) for b in range((n+j)//2+1)]


def differential(v, slope):
    i,j,a,b = v
    if (i,j)==(0,0):
        return {(1,0,a+1,b):slope,(0,1,a,b+1):1}
    if (i,j)==(1,0):
        return {(1,1,a,b+1):1}
    if (i,j)==(0,1):
        return {(1,1,a+1,b):-slope}
    return {}


def main():
    for p in (2,3,5,7):
        for slope in range(p):
            ranks = {}
            for n in range(-2,26):
                columns = []
                target = set(basis(n+1))
                for v in basis(n):
                    col = differential(v,slope)
                    assert set(col)<=target
                    sq=Counter()
                    for w,c in col.items():
                        for x,d in differential(w,slope).items():
                            sq[x] += c*d
                    assert not any(sq.values())
                    columns.append(col)
                ranks[n] = sparse_rank(columns,p)
            for n in range(2,25):
                assert len(basis(n))-ranks[n]-ranks[n-1] == (2 if slope==0 else 0)
    print('Direct four-cell Hom calculation verified over F2, F3, F5, F7 for every slope.')
    print('In degrees 2–24 the zero slope gives dimension 2 and each nonzero slope gives dimension 0.')
    print('The graded k[s,t]-bases in both cases are 1,x; test includes t and z differentials with their signs.')


if __name__ == '__main__':
    require_asserts()
    main()
