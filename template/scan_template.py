"""Arithmetic scan of the half-Steiner template.

Data: tiles T(m) decomposing K_{P x (m-2)}; adjoining a point gives S(2, m-1, 1+P(m-2)).
Special point: Hadamard order P+1.  Ordinary points: Hadamard order r+1 = (P+1)/2, r = (P-1)/2
tiles per point.  Frame: N = (P+1)(P(m-2)+2)/2 unit-norm-P columns, rank P + t(m-1),
t = P(P-1)(m-2)/(m(m-1)) tiles; Seidel eigenvalues rho1 = (P+1)m/2 - P and rho2 = -P.
A sign vector in an eigenspace gives a regular graph of degree (N-1-rho)/2, hence an
srg(N, k, lam, mu) and its complement.
"""
import csv
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
TABLE = HERE.parent / 'data' / 'brouwer_table_snapshot.csv'

def hadamard(n):
    return n in (1, 2) or n % 4 == 0

def steiner_admissible(k, V):
    return (V - 1) % (k - 1) == 0 and (V * (V - 1)) % (k * (k - 1)) == 0

def srg_from_regular(N, rho_on_ones, rho1, rho2):
    """Regular graph in a regular two-graph whose all-ones Seidel eigenvalue is rho_on_ones."""
    k = (N - 1 - rho_on_ones) // 2
    r, s = (-1 - rho2) // 2, (-1 - rho1) // 2   # adjacency restricted eigenvalues
    lam = k + r + s + r * s
    mu = k + r * s
    return (N, k, lam, mu)

def main(NMAX=1300):
    table = {}
    with TABLE.open() as f:
        for row in csv.DictReader(f):
            key = tuple(int(row[c]) for c in ('v', 'k', 'lam', 'mu'))
            table[key] = (row['status'], row['comments'])
    print(f'{"P":>3} {"m":>3} {"S(2,k,V)":>12} {"N":>6} {"rank":>5} {"rho1":>5} {"rho2":>5}  members')
    for P in range(3, 200, 2):
        if not (hadamard(P + 1) and hadamard((P + 1) // 2)):
            continue
        for m in range(3, 40):
            k, V = m - 1, 1 + P * (m - 2)
            N = (P + 1) * (P * (m - 2) + 2) // 2
            if N > NMAX:
                break
            tiles = Fraction(P * (P - 1) * (m - 2), m * (m - 1))
            if tiles.denominator != 1 or not steiner_admissible(k, V):
                continue
            rank = P + int(tiles) * (m - 1)
            rho1, rho2 = (P + 1) * m // 2 - P, -P
            assert rho1 * rho2 == -(N - 1) and rank * rho1 + (N - rank) * rho2 == 0
            out = []
            for rho in (rho1, rho2):
                for prm in (srg_from_regular(N, rho, rho1, rho2),):
                    comp = (N, N - 1 - prm[1], N - 2 - 2 * prm[1] + prm[3], N - 2 * prm[1] + prm[2])
                    st = table.get(prm) or table.get(comp) or ('not listed', '')
                    out.append(f'{prm}[{st[0]}] {st[1][:40]}')
            print(f'{P:>3} {m:>3} {f"S(2,{k},{V})":>12} {N:>6} {rank:>5} {rho1:>5} {rho2:>5}  ' + ' | '.join(out))

if __name__ == '__main__':
    main()
