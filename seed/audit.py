"""Independent standard-library audit of the files and the 15-entry certificate."""
from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path

HERE = Path(__file__).resolve().parent

def read_graph(path):
    data=path.read_bytes().strip()
    assert data[0]==126 and all(63<=k<=126 for k in data)
    n=0
    for c in data[1:4]:
        n=n*64+c-63
    assert n==496 and len(data)==4+(n*(n-1)//2+5)//6
    A=[0]*n
    pairs=((i,j) for j in range(n) for i in range(j))
    for c in data[4:]:
        for shift in (5,4,3,2,1,0):
            i,j=next(pairs)
            if (c-63) & (1<<shift):
                A[i] |= 1<<j
                A[j] |= 1<<i
    return A

def verify_graph(A,k,l,m):
    assert len(A)==496 and all(r.bit_count()==k and not(r>>i&1) for i,r in enumerate(A))
    counts=Counter((A[i]>>j&1,(A[i]&A[j]).bit_count()) for i,j in combinations(range(len(A)),2))
    edges=496*k//2
    assert counts=={(1,l):edges,(0,m):122760-edges}
    return sorted(counts.items())

def independent_witness():
    rows=[tuple(map(int,q.split())) for q in (HERE/'starter.tsv').read_text().splitlines() if q and q[0]!='#']
    vertex_flags=[[] for _ in range(60)]
    design={frozenset([60,*range(i,i+4)]) for i in range(0,60,4)}
    for i,j in product(range(7),range(4)):
        t=4*i+j
        stars=[set() for _ in range(6)]
        for a,b,s,v in rows:
            y=(v%4)^j
            p=56+y if s==2 else 28*s+4*((v//4+i)%7)+y
            stars[a].add(p)
            stars[b].add(p)
            vertex_flags[p].append((t,a,b))
        assert all(len(l)==5 for l in stars)
        design.update(map(frozenset,stars))
    exported=[frozenset(map(int,s.split())) for s in (HERE/'steiner_2_61_5.txt').read_text().splitlines()]
    assert len(exported)==len(set(exported))==183 and set(exported)==design
    paircount=Counter(p for b in design for p in combinations(sorted(b),2))
    assert paircount==Counter(combinations(range(61),2))
    # Direct signed inner products of sparse columns, including their switching signs.
    sign=lambda a,b:(-1)**((a&b).bit_count() % 2)
    V=[]
    for x in range(16):
        V.append({('g',j):sign(j+1,x)*sign(15,x) for j in range(15)})
    for p,flags in enumerate(vertex_flags):
        flags.sort(key=lambda v:((v[1]<2)!=(v[2]<2),v[0]))
        assert len(flags)==7
        for x in range(8):
            w=1 if p>=56 else (sign(4,x)+sign(5,x)+sign(6,x)-sign(7,x))//2
            assert abs(w)==1
            v={('g',p//4):w}
            for c,(t,a,b) in enumerate(flags,1):
                value=w*sign(c,x)
                if c==7 and p<56:value=-value
                v[t,a]=value
                v[t,b]=-value
            V.append(v)
    assert len(V)==496 and all(sum(z*z for z in v.values())==15 for v in V)
    A=[0]*len(V)
    for i,j in combinations(range(len(V)),2):
        ip=sum(a*V[j].get(c,0) for c,a in V[i].items())
        assert ip in (1,-1)
        if ip==-1:
            A[i] |= 1<<j
            A[j] |= 1<<i
    return A

def main():
    graphs=[]
    for k,l,m in [(231,102,112),(264,144,136)]:
        path=HERE/f'srg_496_{k}_{l}_{m}.g6'
        A=read_graph(path)
        print(path.name, sha256(path.read_bytes()).hexdigest())
        print('all vertex degrees',k,'all unordered pairs',verify_graph(A,k,l,m))
        graphs.append(A)
    mask=(1<<496)-1
    assert all(a^b==mask^(1<<i) for i,(a,b) in enumerate(zip(*graphs)))
    assert independent_witness()==graphs[0]
    print('Starter regenerated every edge. The second file is the exact complement.')
    print('Steiner 2-(61,5,1) verified on all 183 blocks and all 1830 point pairs.')
if __name__=='__main__':main()
