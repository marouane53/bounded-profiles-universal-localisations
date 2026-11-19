"""Shared exact finite-field and word utilities for the support-theory checks."""
from collections import Counter
from itertools import product
from functools import cache
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from verify_profile_classification import sparse_rank  # exact sparse column elimination


def rank(rows, prime):
    return sparse_rank([{j: x for j, x in enumerate(row) if x} for row in rows], prime)


def tidy(poly, prime=None):
    if prime:
        return {w: c % prime for w, c in poly.items() if c % prime}
    return {w: c for w, c in poly.items() if c}


def scaled_add(out, vec, c=1):
    for key, value in vec.items():
        out[key] += c * value


def sign(n):
    return -1 if n % 2 else 1


@cache
def weighted_words(weights, d):
    if d < 0:
        return ()
    if not d:
        return ((),)
    return tuple((i,) + tail for i, step in enumerate(weights)
                 for tail in weighted_words(weights, d-step))


def mul(a, b, p):
    return tuple(tuple(sum(x*y for x,y in zip(row,col)) % p for col in zip(*b)) for row in a)


def ident(n):
    return tuple(tuple(int(i == j) for j in range(n)) for i in range(n))


def qpower(q, e, p):
    out = ident(len(q))
    while e:
        if e % 2:
            out = mul(out, q, p)
        q = mul(q, q, p)
        e //= 2
    return out


def require_asserts():
    if not __debug__:
        raise SystemExit('Run with Python assertions enabled.')
