"""Check signed double cuts on finite dg cells, and the recovery map on normal words."""
from collections import Counter
from itertools import product, combinations
from random import Random
from check_lib import scaled_add, sign, tidy, require_asserts


def deriv(w):
    # All letters have degree one; d(t)=xx and d of every other letter is zero.
    out = Counter()
    for i, c in enumerate(w):
        if c == 't':
            out[w[:i]+'xx'+w[i+1:]] += sign(i)
    return tidy(out)


def tensor_d(i, w, j):
    out = Counter({(i,v,j): sign(i)*c for v,c in deriv(w).items()})
    if i:
        out[0,'zz'+w,j] += 1
    if j:
        out[i,w+'zz',0] += sign(i+len(w))
    return tidy(out)


def cuts(i, w, j):
    if 'b' in w:             # Quotient kills b and leaves d(t)=xx.
        return {}
    return {(i,w[:a],w[a+1:b],w[b+1:],j):1
            for a,b in combinations([n for n,c in enumerate(w) if c=='z'],2)}


def target_d(i, u, v, w, j):
    out = Counter()
    if i:
        out[0,'zz'+u,v,w,j] += 1
    for a,c in deriv(u).items():
        out[i,a,v,w,j] += sign(i)*c
    for a,c in deriv(v).items():
        out[i,u,a,w,j] += sign(i+len(u)+1)*c
    for a,c in deriv(w).items():
        out[i,u,v,a,j] += sign(i+len(u)+len(v)+2)*c
    if j:
        out[i,u,v,w+'zz',0] += sign(i+len(u)+len(v)+len(w)+2)
    return tidy(out)


def signed_check():
    rng = Random(142)
    words = [''.join(w) for n in range(6) for w in product('xztb', repeat=n)]
    words += [''.join(rng.choice('xzztt') for _ in range(n)) for n in range(6,13) for _ in range(100)]
    nonzero_before = 0
    for w in words:
        for i,j in product(range(2),repeat=2):
            error = Counter()
            for e,c in cuts(i,w,j).items():
                scaled_add(error, target_d(*e), c)
            for e,c in tensor_d(i,w,j).items():
                scaled_add(error, cuts(*e), -c)
            error = tidy(error)
            nonzero_before += bool(error)
            # This checks a stronger assertion before quotienting cutoff boundaries.
            assert all(a+len(u)<2 or len(y)+b<2 for a,u,v,y,b in error)
    assert nonzero_before
    print(f'Over Z: {4*len(words)} signed chain identities after outer cutoffs, '
          f'{nonzero_before} had actual endpoint errors before truncation.')


def recovery():
    total = 0
    for n in (2,3,4):
        def normal(d):
            return [''.join(w) for w in product('xz',repeat=d) if 'x'*n not in ''.join(w)]
        for p,q in product((2,3),repeat=2):
            for m in range(5):
                for u,v,w in product([s for s in normal(p) if not s.startswith('zz')], normal(m),
                                     [s for s in normal(q) if not s.endswith('zz')]):
                    inserted = u+'z'+v+'z'+w
                    assert 'x'*n not in inserted and not inserted.startswith('zz') and not inserted.endswith('zz')
                    kept = {(i,a,b,c,j):k for (i,a,b,c,j),k in cuts(0,inserted,0).items()
                            if len(a)>=p and len(b)>=m and len(c)>=q}
                    assert kept == {(0,u,v,w,0):1}
                    total += 1
    print(f'{total} outer-class/middle-word triples recovered in fibres x^n=0, n=2,3,4; middle degrees 0–4.')


if __name__ == '__main__':
    require_asserts()
    signed_check()
    recovery()
