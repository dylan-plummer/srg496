# Explicit graphs with the order-84 symmetry

## The nine seeds

Write the ordinary long sheets as \(C_7\times\mathbb F_2^2\), as in the
original certificate. Use \(\mathbb F_8=\mathbb F_2[\alpha]/(\alpha^3+
\alpha+1)\) for the eight coordinates within an ordinary cell, and use
\(\chi_c(x)=(-1)^{\operatorname{Tr}(cx)}\). On the two long sheets,
at the origin, label each of the seven incident tile flags by a
different nonzero \(c\). Impose two conditions.

* The crossing flags have \(\operatorname{Tr}(c)=1\).
* The automorphism doubling the 7-coordinate and cycling the nonzero
  binary colors sends a flag labeled \(c\) to a flag labeled \(c^2\).

There are nine choices on each sheet. Transport them by translation.
The flags of each short-sheet point use \(c=\alpha^i\) on the tile
whose 7-coordinate is \(i\). For the special cell write its arguments
as \((x,e)\in\mathbb F_8\times\mathbb F_2\). The group coordinates
on long sheet \(s\) and position \(i\) are
\(\chi_{\alpha^i}(x)(-1)^{es}\), and the short group coordinate is
\((-1)^e\). Reverse the long crossing flag with \(c=1\). Orient the
triangle on moving symbols cyclically. See [family.py](family.py).

Translations preserve the raw Seidel matrix. A translation by \((r,b)\)
leaves the inner coordinate on long cells alone and multiplies the inner
field coordinate on the short and special cells by \(\alpha^{-r}\).
The order-three automorphism acts by Frobenius on the inner field
coordinates and by the cycle \((3,4,5)\) on tile symbols. The resulting
group is \((C_7\times C_2^2)\rtimes C_3\), with the indicated nontrivial
actions on both factors. Four explicit generators preserving the raw
matrix and the regular switch are checked in
[balanced_census.py](balanced_census.py).

The 81 choices give nine isomorphism classes. For each, the uncolored
degree-231 graph and its regular two-graph both have full group order
84. To obtain the latter order, partition the vertices by the canonical
isomorphism class of the descendant at that vertex. Descendants are
isomorphic exactly when their roots are conjugate in the two-graph.
The automorphism group of a descendant is the stabilizer of its root,
since the isolated root can be reinserted uniquely. Orbit size times
stabilizer order gives 84 for every root class. Testing a representative
of each known graph-group orbit covers all roots. The nine sets of
descendants are disjoint, giving 14 classes per seed and 126 in total.
Bliss computed the first table and nauty checked the exported graphs.

## Switching in the positive eigenspace

Pick one group coordinate to have value one. In each six-symbol tile
choose a subset and put signed one-half on its symbols, zero on the
others. Normalize the subset to omit symbol zero, using the freedom
to add a constant to tile potentials. At a point outside the chosen
group the set of crossing flags must be four characters of an affine
plane avoiding zero. Thus the crossing bits must be the values of
a nonzero binary linear functional on the seven characters. They
vanish at points in the chosen group. These are linear equations for
the tile cuts, followed by a nonzero test. For each potential cut
solution, the product −1 condition on the four crossing signs gives
a second binary linear system for the tile signs. Zero tile cuts
are assigned a fixed sign.

The raw word is \(F^Tz\) and uniquely determines these normalized
potential data. Permuting the words by the full group of order 84
therefore classifies these switches. There are 64 words on eight
seeds and 76 on one. The totals by full graph automorphism order are

| Order | Graph classes |
|---:|---:|
| 2 | 9 |
| 12 | 18 |
| 14 | 11 |
| 84 | 18 |

All 56 graphs are individually certified.
[Compressed graph6 stream](graphs/all_496_231.g6.gz) ·
[cut representatives](cut_classes.tsv.gz) ·
[canonical certificates](all_496_231_audit.tsv).
The 126 descendants are in a second
[graph6 stream](graphs/all_495_238.g6.gz) in the order of
[the descendant table](descendants.tsv).

## All kernel words constant on each ordinary cell

If the special cell carries the constant +1 word, two of four ordinary
cells in each group must carry −1. There are \(6^{15}\) labeled choices.
For a general special sign word the kernel equations say its nonconstant
Walsh coefficients equal −8 times the sums of the four ordinary signs
in the corresponding groups. A nonzero such coefficient must be ±16,
so its absolute value is maximal and the special word is a single
character up to sign. Normalize that sign to be positive. Then one
group has three negative cells and the other fourteen each have two.
This gives \(15\cdot4\cdot6^{14}\) further words. The normalization
already identifies antipodal words; no two normalized words give the
same labeled switching.

The order-84 group acts on these words by ordinary permutation, since
it preserves the raw Seidel matrix. Conversely any isomorphism of
switches must belong to the two-graph group. Thus Burnside's lemma
counts the unlabeled graphs exactly within this family. The two
fixed-word counts of an element, and their frequencies, are

| Constant star | Character star | Elements |
|---:|---:|---:|
| 470184984576 | 4701849845760 | 1 |
| 32768 | 0 | 3 |
| 216 | 144 | 6 |
| 8 | 0 | 18 |
| 0 | 0 | 56 |

There are 61,571,844,416 degree-240 graphs per seed and 554,146,599,744
across the nine. This is an exact count in a restricted family of
switches on nine specified two-graphs. The larger lower bound in
[the proof](PROOF.md) permits all the flag decorations and three
supporting designs.

## Reproduction

Run the following in the srg environment from this directory.

```sh
python census.py
python companions.py
python switch_classes.py
python balanced_census.py
python export.py
python audit_exports.py --canonical
g++ -O3 -std=c++20 starters/enumerate.cpp -o /tmp/srg496_starters
/tmp/srg496_starters > /tmp/srg496_starters.txt
python starters/classify.py /tmp/srg496_starters.txt
python starters/recover.py
python starters/blind_audit.py
python starters/plane_count.py
```

The export audit decodes graph6 directly and checks the degree and
neighborhood intersection of every unordered pair. The optional
canonical step calls nauty to verify pairwise nonisomorphism and the
reported full groups. The seed constructor reads the fifteen entries in
the original sibling directory. The two alternate starter files and
their block lists are included in this directory.
