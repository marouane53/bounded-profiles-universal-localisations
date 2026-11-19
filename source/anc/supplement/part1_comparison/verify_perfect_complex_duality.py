#!/usr/bin/env python3
"""Exact sign check for the chain-level perfect-complex duality.

Over GF(101), construct nontrivial three-term cochain complexes X and Y and
verify degree by degree that

    Y tensor Hom(X,k)  -->  Hom(X,Y),  y tensor f |-> (x |-> f(x)y)

is an isomorphism of complexes with the cohomological Hom and tensor
differentials used in the manuscript.  A deliberately wrong Koszul twist is
also tested and must fail.  No floating point arithmetic is used.
"""
from __future__ import annotations

P = 101


def mm(A, B):
    if not A or not B:
        return []
    return [[sum(A[i][t] * B[t][j] for t in range(len(B))) % P
             for j in range(len(B[0]))] for i in range(len(A))]


def transpose(A):
    return [list(row) for row in zip(*A)] if A else []


def zeros(r, c):
    return [[0] * c for _ in range(r)]


def block_matrix(row_blocks):
    heights = [len(row[0]) if row and row[0] else 0 for row in row_blocks]
    out = []
    for rb, h in zip(row_blocks, heights):
        for i in range(h):
            line = []
            for block in rb:
                line.extend(block[i] if block else [])
            out.append(line)
    return out


def rank(A):
    A = [row[:] for row in A]
    if not A:
        return 0
    m, n = len(A), len(A[0])
    r = 0
    for c in range(n):
        q = next((i for i in range(r, m) if A[i][c] % P), None)
        if q is None:
            continue
        A[r], A[q] = A[q], A[r]
        inv = pow(A[r][c] % P, P - 2, P)
        A[r] = [(inv * x) % P for x in A[r]]
        for i in range(m):
            if i != r and A[i][c] % P:
                a = A[i][c] % P
                A[i] = [(A[i][j] - a * A[r][j]) % P for j in range(n)]
        r += 1
        if r == m:
            break
    return r


# X^{-1}->X^0->X^1 and Y^{-1}->Y^0->Y^1.
Xdims = {-1: 2, 0: 3, 1: 2}
Ydims = {-1: 1, 0: 2, 1: 1}
dX = {
    -1: [[1, 0], [0, 0], [0, 0]],
     0: [[0, 1, 0], [0, 0, 1]],
}
dY = {
    -1: [[1], [0]],
     0: [[0, 1]],
}
assert mm(dX[0], dX[-1]) == zeros(2, 2)
assert mm(dY[0], dY[-1]) == zeros(1, 1)


def hom_basis(n):
    # tuple (source degree p, target degree p+n, target coordinate, source coordinate)
    ans = []
    for p, dx in sorted(Xdims.items()):
        q = p + n
        if q in Ydims:
            for i in range(Ydims[q]):
                for j in range(dx):
                    ans.append((p, q, i, j))
    return ans


def tensor_basis(n):
    # y has degree a, f has degree b and evaluates X^{-b}; total a+b=n.
    ans = []
    for a, dy in sorted(Ydims.items()):
        for b in range(-1, 2):
            p = -b
            if p in Xdims and a + b == n:
                for i in range(dy):
                    for j in range(Xdims[p]):
                        ans.append((a, b, i, j))
    return ans


def phi(n, wrong=False):
    tb = tensor_basis(n); hb = hom_basis(n)
    hi = {x: i for i, x in enumerate(hb)}
    M = zeros(len(hb), len(tb))
    for c, (a, b, i, j) in enumerate(tb):
        key = (-b, a, i, j)
        sign = -1 if wrong and (a * b) % 2 else 1
        M[hi[key]][c] = sign % P
    return M


def hom_differential(n):
    src, tgt = hom_basis(n), hom_basis(n + 1)
    ti = {x: i for i, x in enumerate(tgt)}
    D = zeros(len(tgt), len(src))
    for c, (p, q, i, j) in enumerate(src):
        # d_Y o g
        if q in dY:
            A = dY[q]
            for ip in range(len(A)):
                if A[ip][i] % P:
                    D[ti[(p, q + 1, ip, j)]][c] = (D[ti[(p, q + 1, ip, j)]][c] + A[ip][i]) % P
        # -(-1)^n g o d_X
        if p - 1 in dX:
            A = dX[p - 1]
            coeff = (-1 if n % 2 == 0 else 1) % P
            for jp in range(len(A[0])):
                if A[j][jp] % P:
                    D[ti[(p - 1, q, i, jp)]][c] = (D[ti[(p - 1, q, i, jp)]][c] + coeff * A[j][jp]) % P
    return D


def tensor_differential(n):
    src, tgt = tensor_basis(n), tensor_basis(n + 1)
    ti = {x: i for i, x in enumerate(tgt)}
    D = zeros(len(tgt), len(src))
    for c, (a, b, i, j) in enumerate(src):
        # d_Y y tensor f
        if a in dY:
            A = dY[a]
            for ip in range(len(A)):
                if A[ip][i] % P:
                    D[ti[(a + 1, b, ip, j)]][c] = (D[ti[(a + 1, b, ip, j)]][c] + A[ip][i]) % P
        # (-1)^a y tensor d(f), and d(f)=(-1)^{b+1} f d_X.
        p = -b
        if p - 1 in dX:
            A = dX[p - 1]
            coeff = (-1 if (a + b + 1) % 2 else 1) % P
            for jp in range(len(A[0])):
                if A[j][jp] % P:
                    D[ti[(a, b + 1, i, jp)]][c] = (D[ti[(a, b + 1, i, jp)]][c] + coeff * A[j][jp]) % P
    return D


print('Perfect-complex duality certificate over GF(101)')
all_ok = True
wrong_fails = False
for n in range(-2, 3):
    Pn = phi(n)
    assert len(Pn) == len(hom_basis(n)) and (not Pn or len(Pn[0]) == len(tensor_basis(n)))
    assert rank(Pn) == len(hom_basis(n)) == len(tensor_basis(n))
    if n < 2:
        lhs = mm(phi(n + 1), tensor_differential(n))
        rhs = mm(hom_differential(n), phi(n))
        ok = lhs == rhs
        all_ok &= ok
        wrong_lhs = mm(phi(n + 1, wrong=True), tensor_differential(n))
        wrong_rhs = mm(hom_differential(n), phi(n, wrong=True))
        wrong_fails |= wrong_lhs != wrong_rhs
        print(f'degree {n:2d}: dim={len(hom_basis(n)):2d}, chain square={"PASS" if ok else "FAIL"}')
    else:
        print(f'degree {n:2d}: dim={len(hom_basis(n)):2d}, terminal isomorphism=PASS')
assert all_ok
assert wrong_fails
print('deliberate sign twist (-1)^(|y||f|) fails a chain square: PASS')
print('ALL CHECKS PASSED')
