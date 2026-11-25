"""Audit the released graph6 streams using integer neighbourhood intersections."""
from collections import Counter
import csv
import gzip
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import sys


HERE=Path(__file__).resolve().parent


def decode(line):
    data=line.strip()
    assert data[0]==126 and all(63<=c<=126 for c in data)
    n=(data[1]-63)*4096+(data[2]-63)*64+data[3]-63
    assert len(data)==4+(n*(n-1)//2+5)//6
    rows=[0]*n
    bits=((c-63)>>s&1 for c in data[4:] for s in (5,4,3,2,1,0))
    for j in range(n):
        for i in range(j):
            if next(bits):
                rows[i] |= 1<<j
                rows[j] |= 1<<i
    assert not any(bits)
    return rows


def check(A,n,k,lam,mu):
    assert len(A)==n and all(r.bit_count()==k for r in A)
    for i,j in combinations(range(n),2):
        expected=lam if A[i]>>j&1 else mu
        assert (A[i]&A[j]).bit_count()==expected,(i,j,expected)


def nauty(A):
    import pynauty
    G=pynauty.Graph(len(A),adjacency_dict={i:[j for j in range(len(A)) if r>>j&1] for i,r in enumerate(A)})
    cert=pynauty.certificate(G)
    _,s,e,_,_=pynauty.autgrp(G)
    return sha256(cert).hexdigest(),int(s*10**e)


def main():
    canon='--canonical' in sys.argv
    for name,params,index,field in [
        ('all_496_231',(496,231,102,112),'cut_classes.tsv.gz','automorphisms'),
        ('all_495_238',(495,238,109,119),'descendants.tsv','point_stabilizer')]:
        f=gzip.open(HERE/index,'rt') if index.endswith('.gz') else (HERE/index).open()
        with f as manifest:
            records=list(csv.DictReader(manifest,delimiter='\t'))
        path=HERE/'graphs'/f'{name}.g6.gz'
        with gzip.open(path,'rb') as graph_file:
            lines=graph_file.readlines()
        assert len(lines)==len(records)
        seen=set()
        sizes=Counter()
        with (HERE/f'{name}_audit.tsv').open('w') as out:
            print('index\tgraph6_sha256\tnauty_sha256\tautomorphisms',file=out)
            for i,(line,record) in enumerate(zip(lines,records),1):
                A=decode(line)
                check(A,*params)
                cert,order=('-',0)
                if canon:
                    cert,order=nauty(A)
                    assert order==int(record[field])
                    assert cert not in seen
                    seen.add(cert)
                    sizes[order]+=1
                print(i,sha256(line).hexdigest(),cert,order,sep='\t',file=out,flush=True)
        print(name, len(lines),'graphs',len(lines)*params[0]*(params[0]-1)//2,'unordered pairs, group orders',sorted(sizes.items()),flush=True)
    comp=sorted((HERE/'graphs').glob('496_240_*.g6'))
    for path in comp:
        check(decode(path.read_bytes()),496,240,120,112)
    print('Degree 240:',len(comp),'examples, all neighbourhoods checked.')
    for path in sorted((HERE/'graphs').glob('496_231_*.g6')):
        check(decode(path.read_bytes()),496,231,102,112)
    for path in sorted((HERE/'graphs').glob('495_238_*.g6')):
        check(decode(path.read_bytes()),495,238,109,119)
    print('Individual seed and descendant files verified.')


if __name__=='__main__':
    main()
