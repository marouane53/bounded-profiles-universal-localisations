"""Small exact linear-algebra and finite-algebra utilities over Q.
No numerical packages or floating-point operations are used.
"""
from fractions import Fraction as F
from itertools import product


def vec(n, j=None):
    v = [F(0) for _ in range(n)]
    if j is not None: v[j] = F(1)
    return v


def add(x, y): return [a+b for a,b in zip(x,y)]
def sub(x, y): return [a-b for a,b in zip(x,y)]
def scale(a, x): return [F(a)*b for b in x]
def lincomb(cols, coeff, n=None):
    if n is None: n = len(cols[0]) if cols else 0
    out = vec(n)
    for c,v in zip(coeff,cols):
        if c: out = add(out,scale(c,v))
    return out

def transpose(a): return [list(x) for x in zip(*a)] if a else []
def columns(a): return transpose(a)
def hstack(cols, n=None):
    if n is None: n=len(cols[0]) if cols else 0
    return [[v[i] for v in cols] for i in range(n)]
def matvec(a,x): return [sum((F(c)*d for c,d in zip(row,x)),F(0)) for row in a]
def matmul(a,b):
    bt=transpose(b)
    return [[sum((x*y for x,y in zip(row,col)),F(0)) for col in bt] for row in a]

def rref(a,ncols=None):
    a=[[F(x) for x in row] for row in a]
    if ncols is None: ncols=len(a[0]) if a else 0
    piv=[]; r=0
    for c in range(ncols):
        p=next((i for i in range(r,len(a)) if a[i][c]),None)
        if p is None: continue
        a[r],a[p]=a[p],a[r]
        q=a[r][c]; a[r]=[x/q for x in a[r]]
        for i in range(len(a)):
            if i != r and a[i][c]:
                q=a[i][c]; a[i]=[x-q*y for x,y in zip(a[i],a[r])]
        piv.append(c); r+=1
        if r==len(a): break
    return a,piv

def rank(a,ncols=None): return len(rref(a,ncols)[1])
def independent(cols,n=None):
    ans=[]; r=0
    for v in cols:
        nr=rank(hstack(ans+[v],n))
        if nr>r: ans.append(v);r=nr
    return ans

def nullspace(a,ncols=None):
    if ncols is None: ncols=len(a[0]) if a else 0
    rr,piv=rref(a,ncols)
    ans=[]
    for c in range(ncols):
        if c in piv: continue
        v=vec(ncols,c)
        for i,p in enumerate(piv): v[p]=-rr[i][c]
        ans.append(v)
    return ans

def solve(a,b,ncols=None):
    if ncols is None: ncols=len(a[0]) if a else 0
    if len(a)!=len(b): raise ValueError('Incompatible equation dimensions')
    rr,piv=rref([row+[y] for row,y in zip(a,b)],ncols)
    if any(all(not x for x in row[:ncols]) and row[ncols] for row in rr):
        return None
    v=vec(ncols)
    for i,p in enumerate(piv): v[p]=rr[i][ncols]
    return v

def coords(cols,v):
    out=solve(hstack(cols,len(v)),v,len(cols))
    if out is None: raise ValueError('Vector not in indicated span')
    return out


class Algebra:
    def __init__(self,names,products,unit):
        self.names=list(names); self.n=len(names)
        self.table={}
        for (i,j),v in products.items():
            if isinstance(v,dict):
                w=vec(self.n)
                for k,c in v.items():w[k]=F(c)
            else:w=[F(c) for c in v]
            self.table[i,j]=w
        self.unit=[F(x) for x in unit]
    def basis(self,i):return vec(self.n,i)
    def mul(self,x,y):
        out=vec(self.n)
        for i,a in enumerate(x):
            if not a:continue
            for j,b in enumerate(y):
                if b and (i,j) in self.table:
                    out=add(out,scale(a*b,self.table[i,j]))
        return out
    def check(self):
        bs=[self.basis(i) for i in range(self.n)]
        for x in bs:
            assert self.mul(self.unit,x)==x==self.mul(x,self.unit)
        for x,y,z in product(bs,repeat=3):
            assert self.mul(self.mul(x,y),z)==self.mul(x,self.mul(y,z))
    def left(self,x):return hstack([self.mul(x,self.basis(i)) for i in range(self.n)])
    def right(self,x):return hstack([self.mul(self.basis(i),x) for i in range(self.n)])


class DGA(Algebra):
    def __init__(self,names,degrees,products,unit,d):
        super().__init__(names,products,unit)
        self.degrees=list(degrees);self.d=[[F(x) for x in row] for row in d]
    def check_dga(self):
        self.check()
        assert len(self.d)==self.n and all(len(row)==self.n for row in self.d)
        assert all(not self.d[i][j] or self.degrees[i]==self.degrees[j]+1 for i in range(self.n) for j in range(self.n))
        assert all(not c or self.degrees[i]==0 for i,c in enumerate(self.unit))
        assert all(not x for row in matmul(self.d,self.d) for x in row)
        bs=[self.basis(i) for i in range(self.n)]
        for i,x in enumerate(bs):
            for j,y in enumerate(bs):
                xy=self.mul(x,y)
                assert all(not c or self.degrees[k]==self.degrees[i]+self.degrees[j]
                           for k,c in enumerate(xy))
                lhs=matvec(self.d,xy)
                rhs=add(self.mul(matvec(self.d,x),y),scale((-1 if self.degrees[i] % 2 else 1),self.mul(x,matvec(self.d,y))))
                assert lhs==rhs,(i,j,lhs,rhs)
