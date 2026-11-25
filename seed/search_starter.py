"""Find a T(6) starter with the order-three symmetry on symbols 3,4,5."""
from itertools import combinations

def rotate(edge):
    permutation = (0, 1, 2, 4, 5, 3)
    return tuple(sorted(permutation[i] for i in edge))

def transform(x):
    """Multiply the Z/7 coordinate by two and cycle the nonzero binary coordinates."""
    return 4 * (2 * (x // 4) % 7) + (0, 2, 3, 1)[x % 4]

def difference(x, y):
    return 4 * ((x // 4 - y // 4) % 7) + ((x ^ y) & 3)

def edge_class(first, second):
    sheet, pos = first
    other, dest = second
    if sheet > other:
        return edge_class(second, first)
    if sheet == other:
        if pos // 4 == dest // 4:
            return None
        return (sheet, other, min(difference(pos, dest), difference(dest, pos)))
    if other == 2:
        return (sheet, other, (pos ^ dest) & 3)
    return (sheet, other, difference(dest, pos))

def extend(tile, occupied, orbits):
    if not orbits:
        return tile if len(occupied) == 60 else None
    edges, sheet = orbits[0]
    for pos in range(1, 28):
        trial = dict(tile)
        classes = set(occupied)
        current = pos
        for edge in edges:
            point = (sheet, current)
            if point in trial.values():
                break
            added = [edge_class(point, q) for f,q in trial.items() if set(edge).intersection(f)]
            if None in added or len(set(added)) != len(added) or classes.intersection(added):
                break
            classes.update(added)
            trial[edge] = point
            current = transform(current)
        else:
            hit = extend(trial, classes, orbits[1:])
            if hit:
                return hit
    return None

def search():
    seen = {(0,1)}
    orbits = []
    for b in range(6):
        for a in range(b):
            e = (a,b)
            if e in seen:
                continue
            orbit = [e]
            while rotate(orbit[-1]) != e:
                orbit.append(rotate(orbit[-1]))
            seen.update(orbit)
            orbits.append(orbit)
    for mask in range(1 << len(orbits)):
        colours = {(0,1):2}
        for i,orbit in enumerate(orbits):
            for edge in orbit:
                colours[edge] = (mask >> i) & 1
        if list(colours.values()).count(0) != 7:
            continue
        adjacent = [(colours[e], colours[f]) for e,f in combinations(colours,2) if set(e).intersection(f)]
        if adjacent.count((0,0)) != 12 or sum(set(pair)=={0,2} for pair in adjacent)!=4:
            continue
        fixed = {(0,1):(2,0)}
        free = []
        for orbit in orbits:
            sheet = colours[orbit[0]]
            if len(orbit)==1:
                fixed[orbit[0]] = (sheet,0)
            else:
                free.append((orbit,sheet))
        initial = [edge_class(fixed[e],fixed[f]) for e,f in combinations(fixed,2) if set(e).intersection(f)]
        if None in initial or len(set(fixed.values()))!=len(fixed) or len(set(initial))!=len(initial):
            continue
        hit = extend(fixed,set(initial),free)
        if hit:
            return hit
    raise ValueError('No starter found')

if __name__=='__main__':
    for edge, point in search().items():
        print(*edge,*point)
