"""Recover each Steiner design from a relabelled, switched graph6 input."""
from collections import defaultdict
from itertools import combinations
from pathlib import Path
import random
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from audit_exports import decode,check
import pynauty


class Partition:
    def __init__(self,n):
        self.parent=list(range(n))

    def find(self,a):
        while a!=self.parent[a]:
            self.parent[a]=self.parent[self.parent[a]]
            a=self.parent[a]
        return a

    def join(self,a,b):
        self.parent[self.find(b)]=self.find(a)

    def parts(self):
        bags=defaultdict(list)
        for i in range(len(self.parent)):
            bags[self.find(i)].append(i)
        return list(bags.values())


def project(row,mask):
    r=row&mask
    return min(r,r^mask)


def recover(A):
    n=len(A)
    full=(1<<n)-1
    M=[r | (1<<i) for i,r in enumerate(A)]
    buckets=defaultdict(list)
    for i,j in combinations(range(n),2):
        buckets[project(M[i]^M[j],full)].append((i,j))
    pieces=Partition(n)
    quad_count=0
    for edges in buckets.values():
        for (a,b),(c,d) in combinations(edges,2):
            if len({a,b,c,d})==4:
                pieces.join(a,b)
                pieces.join(a,c)
                pieces.join(a,d)
                quad_count+=1
    assert quad_count//3==980
    assert sorted(map(len,pieces.parts()))==[8]*60+[16]
    special,=[c for c in pieces.parts() if len(c)==16]
    mask=sum(1<<i for i in special)
    parts=defaultdict(list)
    for i in range(n):
        if not mask>>i&1:
            parts[project(A[i],mask)].append(i)
    assert len(parts)==15
    cells=[]
    blocks=[]
    for vertices in parts.values():
        assert len(vertices)==32
        m=sum(1<<i for i in vertices)
        fibres=Partition(len(vertices))
        for a,b in combinations(range(32),2):
            i,j=vertices[a],vertices[b]
            others=m^(1<<i)^(1<<j)
            parity=(A[i]^A[j])&others
            if (A[i]>>j&1 and parity==0) or (not A[i]>>j&1 and parity==others):
                fibres.join(a,b)
        fibre_list=[[vertices[a] for a in ids] for ids in fibres.parts()]
        assert len(fibre_list)==4 and all(len(c)==8 for c in fibre_list)
        blocks.append([60,*range(len(cells),len(cells)+4)])
        cells.extend(fibre_list)
    masks=[sum(1<<a for a in c) for c in cells]
    flags=Partition(60*7)
    labels=[]
    for p in range(60):
        directions={}
        out={}
        for q in range(60):
            if p//4==q//4:
                continue
            vec=project(A[cells[q][0]],masks[p])
            assert all(project(A[v],masks[p])==vec for v in cells[q])
            out[q]=7*p+directions.setdefault(vec,len(directions))
        assert len(directions)==7
        labels.append(out)
    for p,q in combinations(range(60),2):
        if q in labels[p]:
            flags.join(labels[p][q],labels[q][p])
    tiles=flags.parts()
    assert len(tiles)==28
    for tile in tiles:
        pts=[i//7 for i in tile]
        assert len(pts)==len(set(pts))==15
        active=set(tile)
        nb={p:{q for q,r in labels[p].items() if r in active} for p in pts}
        # In T(6) every 5-clique consists of the 5 edges at a single symbol.
        wanted=[c for c in combinations(pts,5) if all(b in nb[a] for a,b in combinations(c,2))]
        assert len(wanted)==6
        blocks.extend(wanted)
    return quad_count//3,blocks


def cert(blocks):
    adj={i:[] for i in range(244)}
    for i,b in enumerate(blocks,61):
        for v in b:
            adj[v].append(i)
            adj[i].append(v)
    return pynauty.certificate(pynauty.Graph(244,adjacency_dict=adj))


def main():
    rng=random.Random(341155)
    for path in sorted(HERE.glob('design*.g6')):
        orig=decode(path.read_bytes())
        check(orig,496,231,102,112)
        ix=list(range(496))
        rng.shuffle(ix)
        side=[rng.randrange(2) for i in ix]
        A=[sum(1<<j for j in range(496) if j!=i and (orig[ix[i]]>>ix[j]&1)^side[j]^side[i]) for i in range(496)]
        quads,blocks=recover(A)
        ref=[list(map(int,l.split())) for l in path.with_suffix('.blocks').read_text().splitlines()]
        assert cert(blocks)==cert(ref)
        print(path.stem,quads,'qualifying quadruples; recovered the design on 61 points,',len(blocks),'blocks',flush=True)


if __name__=='__main__':
    main()
