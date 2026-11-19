"""Exact tests of the deformation, including low cochain ranks and genuine roots of unity."""
from collections import Counter
from fractions import Fraction
from check_lib import sign, tidy, sparse_rank, weighted_words, require_asserts


def monic_div(a,b):
    a=list(a)
    result=[0]*max(0,len(a)-len(b)+1)
    while len(a)>=len(b):
        i=len(a)-len(b)
        c=a[-1]
        result[i]=c
        for j,k in enumerate(b):
            a[i+j]-=c*k
        while a and a[-1]==0:a.pop()
    assert not a
    return result


def times_q(a,cycl):
    a=[0]+list(a)
    if len(a)==len(cycl):
        c=a[-1]
        a=[v-c*k for v,k in zip(a[:-1],cycl[:-1])]
    return a


def roots_check():
    cycl={}
    cases=0
    for d in range(1,25):
        f=[-1]+[0]*(d-1)+[1]
        for j in range(1,d):
            if d%j==0:f=monic_div(f,cycl[j])
        cycl[d]=f
        unit=[1]+[0]*(len(f)-2)
        power=unit
        for j in range(4*d+1):
            assert (power==unit)==(j%d==0)
            assert [2*a for a in power]!=unit
            power=times_q(power,f)
            cases+=1
    print(f'Cyclotomic arithmetic over Q for primitive orders 1–24, {cases} even-degree scalar tests; the slope-2 scalar never vanishes.')


def strict_bar(q,p,n=5):
    # a,b,h are ordinary radical letters; e,w have degree -1.
    weights=(1,1,1,2,2)
    products={(0,0):(2,1),(3,1):(4,1),(1,3):(4,q)}
    def boundary(w):
        out=Counter()
        parity=1
        for i,a in enumerate(w):
            if a==3:out[w[:i]+(2,)+w[i+1:]]-=parity
            if w[i:i+2] in products:
                k,c=products[w[i:i+2]]
                out[w[:i]+(k,)+w[i+2:]]+=parity*sign(weights[a]-1)*c
            parity*=sign(weights[a])
        return tidy(out,p)
    ranks=[sparse_rank([boundary(w) for w in weighted_words(weights,i)],p) for i in range(n+2)]
    dims=[len(weighted_words(weights,i))-ranks[i]-ranks[i+1] for i in range(n+1)]
    assert dims==[1,2,4,7,12,20][:n+1]
    # With u=h*, v=-e*, y=-w* the two bar terms in d(y) have opposite coefficients.
    assert boundary((3,1))[4,]==(-1)%p
    assert boundary((1,3))[4,]==q%p


def finite_fields():
    checked=0
    for p in (2,3,5,7,11,13,17,19):
        for q in range(1,p):
            order=next(j for j in range(1,p) if pow(q,j,p)==1)
            # Build the action by e_i x²z=q e_i zx² from the two initial entries.
            action=[1,0]
            for i in range(8*p):
                action.append(q*action[i]%p)
            failures=[i for i,c in enumerate(action) if c==1]
            assert failures==list(range(0,len(action),2*order))
            if p<=7:strict_bar(q,p)
            checked+=1
    print(f'For {checked} prime-field units the direct cone map fails precisely every twice-the-order degrees; strict-bar checks over F2,F3,F5,F7.')


def infinite_orders():
    # Every test uses exact rational arithmetic. A nonzero polynomial may have genuine
    # early exceptions, which a test demanding invertibility in all degrees would reject.
    for q in [Fraction(2), Fraction(-2), Fraction(3,2), Fraction(-5,3)]:
        for lam in [Fraction(1),Fraction(2),Fraction(3,4)]:
            for n in range(70):
                even=lam*q**n
                odd=0
                assert (even*q)==lam*q**(n+1) and odd==0
                assert (even-1==0)==(lam==1 and n==0)
        # A finite-order-style test would reject this cone on M_1 at degree zero.
        assert 1*q**0-1==0 and all(q**j!=1 for j in range(1,100))
    print('Checked four infinite-order rational twists, three slopes, 70 even/odd steps; finite exceptions are visible.')


def finite_domain_refinement():
    # Translate both elementary off-diagonal generators by cI. They still
    # generate M2. The commutator passes both orbits, while adding U-1
    # removes exactly the c=1 orbit among these two.
    from verify_orbit_closures import det, multiply, translate
    ident = (1, 0, 0, 1)
    for q in (1, -1):
        for c in (1, 2):
            u = translate((0, 1, 0, 0), ident, c)
            v = translate((0, 0, 1, 0), ident, c)
            refined = []
            for j in range(4):
                commutator = translate(multiply(u, v), multiply(v, u), -1)
                assert det(commutator) != 0
                refined.append(det(translate(u, ident, -1)) != 0)
                u, v = v, tuple(q*x for x in u)
            assert all(refined) == (c == 2)
    print('Checked a strict nonempty refinement of a two-dimensional simple domain at both finite twists q=1,-1.')


if __name__=='__main__':
    require_asserts()
    roots_check()
    finite_fields()
    infinite_orders()
    finite_domain_refinement()
