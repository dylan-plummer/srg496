"""Export all three catalogues as newline-delimited graphs in graph6 format."""
import csv
import gzip
from itertools import groupby
import networkx as nx
import numpy as np
from family import HERE, make
from companions import descendant
from switch_classes import graph_from_code


def graph6(A):
    n=len(A)
    header=bytes((126,63+(n>>12),63+((n>>6)&63),63+(n&63)))
    u,v=np.tril_indices(n,-1)
    bits=np.asarray(A[u,v],dtype=np.uint8)
    bits=np.pad(bits,(0,(-len(bits))%6))
    chars=63+bits.reshape(-1,6)@np.array([32,16,8,4,2,1],dtype=np.uint8)
    return header+chars.tobytes()+b'\n'


def main():
    with gzip.open(HERE/'cut_classes.tsv.gz','rt') as inp,gzip.open(HERE/'graphs'/'all_496_231.g6.gz','wb') as out:
        for pal,rows in groupby(csv.DictReader(inp,delimiter='\t'),key=lambda r:r['palette']):
            F=make(*map(int,pal))
            for r in rows:
                code=int(r['group']),bytes.fromhex(r['cut_hex'])
                out.write(graph6(graph_from_code(F,code)))
    with (HERE/'descendants.tsv').open() as inp,gzip.open(HERE/'graphs'/'all_495_238.g6.gz','wb') as out:
        for pal,rows in groupby(csv.DictReader(inp,delimiter='\t'),key=lambda r:r['palette']):
            G=nx.read_graph6(HERE/'graphs'/f'496_231_{pal}.g6')
            A=nx.to_numpy_array(G,dtype=np.int16)
            S=1-np.eye(len(A),dtype=np.int16)-2*A
            for r in rows:
                out.write(graph6(descendant(S,int(r['root']))))


if __name__=='__main__':
    main()
