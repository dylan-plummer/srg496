"""Independent ordinary graph6 checks including blind recovery of the 976-point parents."""
from collections import defaultdict,Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import csv,random
import igraph as ig

HERE=Path(__file__).resolve().parent
PARAMS={'976_480':(976,480,248,224),'975_462':(975,462,205,231),
        '496_240':(496,240,120,112),'496_231':(496,231,102,112),'495_238':(495,238,109,119)}
def decode(data):
    raw=list(x-63 for x in data.strip())
    assert raw[0]==63
    n=4096*raw[1]+64*raw[2]+raw[3]
    assert len(raw)==4+(n*(n-1)//2+5)//6
    rows=[0]*n
    idx=0
    for j in range(n):
        for i in range(j):
            if raw[4+idx//6]>>(5-idx%6)&1:
                rows[i]|=1<<j
                rows[j]|=1<<i
            idx+=1
    return rows

def check(rows,params):
    n,k,l,m=params
    assert len(rows)==n and all(x.bit_count()==k for x in rows)
    for i,a in enumerate(rows):
        for j in range(i):
            assert (a&rows[j]).bit_count()==(l if a>>j&1 else m)
def polys(rows,negative):
    n=len(rows);mask=(1<<n)-1
    # Bits encode negative entries of S-I. E.g. for A=[S>0] they are the
    # bitwise complement of A (including its diagonal).
    words=[(v|(1<<i) if negative else mask^v) for i,v in enumerate(rows)]
    key=lambda x:min(x,mask^x)
    first={};others=defaultdict(list)
    for j,x in enumerate(words):
        for i in range(j):
            k=key(words[i]^x)
            if k in first:
                others[k].append((i,j))
            else:
                first[k]=(i,j)
    comps=[{i} for i in range(n)]
    root=list(range(n))
    deg=Counter()
    four=set()
    for k,other in others.items():
        for p,q in combinations([first[k]]+other,2):
            quad=tuple(sorted(p+q))
            assert len(set(quad))==4
            if quad in four:
                continue
            four.add(quad);deg.update(quad)
            a=root[quad[0]]
            for x in quad[1:]:
                z=root[x]
                if z==a:continue
                comps[a].update(comps[z])
                for j in comps[z]:root[j]=a
                comps[z].clear()
    parts=[sorted(c) for c in comps if c]
    return words,len(four),parts,deg


def recognize(words,parts,deg):
    # This case is used for a distinguished connected 28-plane Hadamard at
    # the special point. Ordinary connected 16-cells have degree >=11.
    special,=[c for c in parts if len(c)==16 and all(deg[i]==7 for i in c)]
    def projection(vertex,cell):
        mask=sum(((words[vertex]>>u)&1)<<j for j,u in enumerate(cell))
        return min(mask,((1<<len(cell))-1)^mask)
    groups=defaultdict(list)
    for c in parts:
        if c==special:continue
        sig={projection(i,special) for i in c}
        assert len(sig)==1
        groups[sig.pop()].append(c)
    assert len(groups)==15
    full=[special]
    for key,cells in groups.items():
        outside=[c for kk,cs in groups.items() if kk!=key for c in cs]
        halves=defaultdict(list)
        for c in cells:
            if len(c)==16:full.append(c);continue
            assert len(c)==8
            halves[tuple(projection(c[0],t) for t in outside)].append(c)
        for twins in halves.values():
            a,b=twins
            full.append(sorted(a+b))
    assert sorted(map(len,full))==[16]*61
    blocks=set()
    for i,c in enumerate(full):
        by=defaultdict(list)
        for j,t in enumerate(full):
            if j!=i:by[projection(t[0],c)].append(j)
        assert len(by)==15 and all(len(v)==4 for v in by.values())
        blocks.update(tuple(sorted([i]+v)) for v in by.values())
    assert len(blocks)==183
    return design_certificate(blocks)

def design_certificate(blocks):
    g=ig.Graph(244,edges=[(p,61+k) for k,b in enumerate(sorted(blocks)) for p in b])
    perm=g.canonical_permutation(color=[0]*61+[1]*183)
    return sha256(repr(sorted(g.permute_vertices(perm).get_edgelist())).encode()).hexdigest()


def norm(rows):
    n=len(rows);mask=(1<<n)-1;a=rows[0]
    return [((r^a^(mask if a>>i&1 else 0))&~(1<<i)) for i,r in enumerate(rows)]

def main():
    with (HERE/'audit.tsv').open('w') as f:
        print('graph\tsha256\tpairs\tplanes\tcomponent_sizes',file=f)
        groups=defaultdict(dict)
        for file in sorted((HERE/'graphs').glob('*.g6')):
            tag,n,k=file.stem.rsplit('_',2);param=n+'_'+k
            data=file.read_bytes();rows=decode(data);check(rows,PARAMS[param])
            groups[tag][param]=file
            count=''
            profile=''
            if param in ('976_480','496_231'):
                words,count,parts,degs=polys(rows,param=='496_231')
                profile=sorted(Counter(map(len,parts)).items())
                if param=='976_480' and ('_h4' in tag or tag=='mixed_all'):
                    cert=recognize(words,parts,degs)
                    d=int(tag[1]) if tag.startswith('d') else 0
                    path=HERE.parent/'starters'/f'design{d:02}.blocks'
                    reference=design_certificate([tuple(sorted(map(int,line.split()))) for line in path.read_text().splitlines()])
                    assert cert==reference
                    # Scramble the labels and apply arbitrary simultaneous row/column signs.
                    rng=random.Random(9621+d)
                    perm=list(range(len(words)))
                    rng.shuffle(perm)
                    flip=rng.getrandbits(len(words));whole=(1<<len(words))-1
                    scrambled=[(sum(((words[u]>>v)&1)<<k for k,v in enumerate(perm))^flip^
                                (whole if flip>>j&1 else 0)) for j,u in enumerate(perm)]
                    # In the positive adjacency orientation, bit rows of S-I are ~adjacency.
                    _,again,comps,degree=polys([whole^x for x in scrambled],False)
                    assert again==count and recognize(scrambled,comps,degree)==cert
                    print(tag,'recovered source design after random switches and relabeling',cert,flush=True)
            print(file.stem,sha256(data).hexdigest(),int(n)*(int(n)-1)//2,count,profile,sep='\t',file=f,flush=True)
        for tag,files in groups.items():
            rows=lambda k:decode(files[k].read_bytes())
            par=rows('976_480')
            sub=rows('496_240')
            assert [r&((1<<496)-1) for r in par[:496]]==sub
            assert [r&((1<<495)-1) for r in rows('975_462')[:495]]==rows('495_238')
            # Seidel sign is reversed for the two degree orientations.
            def opposite(a):
                mask=(1<<len(a))-1
                return [mask^r^(1<<i) for i,r in enumerate(a)]
            assert norm(opposite(sub))==norm(rows('496_231'))
            for a,b in [('976_480','975_462'),('496_240','495_238')]:
                assert [r>>1 for r in norm(opposite(rows(a)))[1:]]==rows(b)
        print('all ordinary pair counts, induced subgraphs, normalized switches and descendants passed for',len(groups),'chains')


if __name__=='__main__':main()
