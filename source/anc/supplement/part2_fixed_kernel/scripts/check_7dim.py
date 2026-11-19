# Exact (rational) check of the 7-dimensional example, pure Python.
# C = k[x,y]/(x^2,y^2) with basis (1, x, y, xy).  M = C/(x), N = C/(y).
from fractions import Fraction as F

def rank(rows):
    A=[[F(v) for v in r] for r in rows]; 
    if not A: return 0
    m,n=len(A),len(A[0]); r=0
    for c in range(n):
        p=next((i for i in range(r,m) if A[i][c]!=0),None)
        if p is None: continue
        A[r],A[p]=A[p],A[r]
        for i in range(m):
            if i!=r and A[i][c]!=0:
                f=A[i][c]/A[r][c]; A[i]=[a-f*b for a,b in zip(A[i],A[r])]
        r+=1
    return r
def mul(A,B): return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def zero(n,m): return [[0]*m for _ in range(n)]

# left multiplication by x and y on C in basis (1,x,y,xy): column j = image of basis vector j
X=[[0,0,0,0],[1,0,0,0],[0,0,0,0],[0,0,1,0]]
Y=[[0,0,0,0],[0,0,0,0],[1,0,0,0],[0,1,0,0]]
assert mul(X,X)==zero(4,4) and mul(Y,Y)==zero(4,4) and mul(X,Y)==mul(Y,X)

# ker(y: C->C) = yC  (so the resolution of N=C/(y) is ... C -y-> C -y-> C -> N, all maps y)
print("rank(y)=%d, dim ker(y)=%d, dim yC=%d  => ker(y)=yC, 1-periodic resolution of N=C/(y)"%(rank(Y),4-rank(Y),rank(Y)))
print("rank(x)=%d, dim ker(x)=%d, dim xC=%d  => ker(x)=xC, 1-periodic resolution of M=C/(x)"%(rank(X),4-rank(X),rank(X)))

# Tor_i^C(M,N): tensor the y-resolution with M=C/(x) (basis 1,y): y acts on M by Ym
Ym=[[0,0],[1,0]]
print("Tor_0^C(M,N) dim =", 2-rank(Ym), "; Tor_i^C(M,N) dim (i>=1) =", (2-rank(Ym))-rank(Ym), "  (expected: k and 0)")
# Tor_i^C(M,M): tensor the x-resolution with M: x acts by 0 on M
Xm=zero(2,2)
print("Tor_0^C(M,M) dim =", 2-rank(Xm), "; Tor_i^C(M,M) dim (i>=1) =", (2-rank(Xm))-rank(Xm), "  (expected: M, M  i.e. dim 2, 2)")

# stable End_C(N), N=C/(y) (basis 1,x):  End_C(N)=N has dim 2.
# Hom_C(N,C) = {c : y c = 0} = ann(y); composite N -> C -> N=C/(y): reduce c mod (y,xy)
annY_basis=[[0,0,1,0],[0,0,0,1]]   # y, xy
assert all(mul(Y,[[v] for v in c])==zero(4,1) for c in annY_basis)
reduced=[c[:2] for c in annY_basis]  # projection C -> C/(y) keeps coordinates of 1 and x
print("images of 1 under N->C->N:", reduced, " => all zero, stable End_C(N) = End_C(N) = k[x]/(x^2), dim 2 (eq. 13)")

# A = [[k, M],[0, C]]: dim 1+2+4 = 7; right module k_A has projective cover e_1A=(k,M) with kernel (0,M)=M_A,
# whose pd over A equals pd_C(M) = infinity since M (dim 2) is not free over the self-injective 4-dim C.
print("dim A =", 1+2+4, "; kernel of e_1A -> k_A is M_A with dim", 2, "not divisible by 4 => pd(k_A)=infinity (bounded extension fails)")
