"""Enumerate order-16 Hadamard matrices by Orrick switches and distinguish pointed designs."""
from collections import Counter,defaultdict
from itertools import combinations
from pathlib import Path
import json
import pynauty

HERE=Path(__file__).resolve().parent
N=16
FULL=(1<<16)-1

def graph_key(adjacency, colours):
    g=pynauty.Graph(len(adjacency),adjacency_dict=adjacency,vertex_coloring=[set(a) for a in colours])
    _,m,e,orbits,_=pynauty.autgrp(g)
    return pynauty.certificate(g), round(m*10**e), orbits

def transpose(h):
    return tuple(sum(((r>>j)&1)<<i for i,r in enumerate(h)) for j in range(N))

def equivalence(h):
    # Rows and their negatives, columns and their negatives.
    a={i:set() for i in range(64)}
    for i,r in enumerate(h):
        for j in range(16):
            s=(r>>j)&1
            for sign in (0,1):
                u=2*i+sign;v=32+2*j+(s^sign)
                a[u].add(v);a[v].add(u)
    return graph_key(a,[range(32),range(32,64)])[0]

def pointed(h):
    a={i:set() for i in range(46)}
    for i,r in enumerate(h[1:]):
        for j in range(16):
            v=16+2*i+((r>>j)&1)
            a[v].add(j);a[j].add(v)
    return graph_key(a,[range(16),range(16,46)])

def normalize(h, row=0):
    x=h[row]
    other=[r^x for i,r in enumerate(h) if i!=row]
    return (0,*[min(r,FULL^r) for r in other])

def verify(h):
    assert len(h)==16 and all((r^s).bit_count()==8 for r,s in combinations(h,2))

def switch_options(h):
    for ids in combinations(range(16),4):
        a,b,c,d=(h[i] for i in ids)
        if a^b^c^d not in (0,FULL):
            continue
        # For a closed quadruple any one of the four fields gives an H-equivalent switch.
        field=FULL&~(a^b)&~(a^c)
        assert field.bit_count()==4
        yield tuple(x^(field if i in ids else 0) for i,x in enumerate(h))

def main():
    walsh=tuple(sum(((i&j).bit_count()%2)<<j for j in range(16)) for i in range(16))
    unpointed=[walsh]
    seen={equivalence(walsh)}
    for h in unpointed:
        for is_col in (False,True):
            t=transpose(h) if is_col else h
            for x in switch_options(t):
                f=transpose(x) if is_col else x
                verify(f)
                sig=equivalence(f)
                if sig not in seen:
                    unpointed.append(f);seen.add(sig)
        print('unpointed found',len(seen),flush=True)
    assert len(unpointed)==5
    seen_p=set()
    examples=[]
    for i,h in enumerate(unpointed):
        for k in range(16):
            x=normalize(h,k)
            sig,order,_=pointed(x)
            if sig in seen_p:
                continue
            seen_p.add(sig)
            tr=transpose(x)
            planes=[c for c in combinations(range(16),4) if not(tr[c[0]]^tr[c[1]]^tr[c[2]]^tr[c[3]])]
            counts=Counter(j for c in planes for j in c)
            comps=[{v} for v in range(16)]
            for c in planes:
                hits=[s for s in comps if set(c)&s]
                comps=[s for s in comps if s not in hits]+[set.union(*hits)]
            data={'pointed':len(examples),'unpointed':i,'row':k,'order':order,'planes':len(planes),
                  'plane_degrees':[counts[j] for j in range(16)],'components':sorted(map(len,comps)),'matrix':list(x)}
            print({k:v for k,v in data.items() if k!='matrix'},flush=True)
            examples.append(data)
    (HERE/'hadamards.json').write_text(json.dumps(examples,indent=1)+'\n')
if __name__=='__main__':main()
