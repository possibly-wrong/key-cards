import numpy as np
from functools import reduce
from collections import Counter
from itertools import chain, combinations

def cycles(p):
    """Convert permutation to list of disjoint cycles."""
    p = dict(enumerate(p))
    while p:
        x = next(iter(p))
        cycle = []
        while x in p:
            cycle.append(x)
            x = p.pop(x)
        yield cycle

def merge(s, p):
    """Merge union-find s with permutation p (as cycles)."""
    def find(x):
        while s[x] != x:
            x = s[x]
        return x
    def union(x, y):
        x = find(x)
        s[find(y)] = x
        return x
    for cycle in p:
        reduce(union, cycle)
    for x in range(len(s)):
        s[x] = find(x)
    return s

def cycle_index_term(s, k=2):
    """Convert union-find s to cycle index monomial at x[i]=k."""
    #return prod(x[i]**j for i, j in Counter(Counter(s).values()).items())
    return k ** sum(Counter(Counter(s).values()).values())

def asymmetric_colorings(group, k=2):
    """Number of k-colorings with no symmetries in the given group."""

    # Group G acts on (colorings of) X = {0, 1, 2, ..., n-1}.
    G = list(group)
    n = sum(len(cycle) for cycle in G[0])

    # Compute inclusion-exclusion sum over subsets of G-e.
    G = [g for g in G if len(g) < n]
    return sum((-1) ** len(subset) *
               cycle_index_term(reduce(merge, subset, list(range(n))), k)
               for subset in chain.from_iterable(combinations(G, r)
                                                 for r in range(len(G) + 1)))

def cycle_index(group, k=2):
    """Cycle index for given group action (as cycles) at x[i]=k."""
    G = list(group)
    n = sum(len(cycle) for cycle in G[0])
    return (sum(cycle_index_term(merge(list(range(n)), g), k) for g in G)
            // len(G))

# Following are examples of group actions.

def dihedral_grid(n):
    """Dihedral group action on squares of nxn grid (as cycles)."""
    assert(n > 1)
    g = np.split(np.arange(n ** 2), n)
    for i in range(2):
        for j in range(4):
            yield list(cycles(np.ndarray.flatten(np.array(g))))
            g = np.rot90(g)
        g = np.transpose(g)

def dihedral_bracelet(n):
    """Dihedral group action on n-bead bracelets (as cycles)."""
    assert(n > 2)
    for i in range(n):
        g = list(range(i, n)) + list(range(i))
        yield list(cycles(g))
        yield list(cycles(reversed(g)))

def cyclic_necklace(n):
    """Cyclic group action on n-bead necklaces (as cycles)."""
    for i in range(n):
        g = list(range(i, n)) + list(range(i))
        yield list(cycles(g))

##import isohedra
##def face_symmetries(vertices, faces):
##    """Action of symmetries of given isohedron on its faces (as cycles)."""
##    G = {g: faces[:] for g in isohedra.symmetry_group(vertices, faces, False)}
##    for i, face_i in enumerate(faces):
##        for j, face_j in enumerate(faces):
##            for R in isohedra.symmetries(vertices, face_i, face_j, False):
##                for g, p in G.items():
##                    if abs(R - g).max() < 1e-3:
##                        p[i] = j
##    return [list(cycles(p)) for p in G.values()]

if __name__ == '__main__':
    for n in range(2, 9):
        print(n,
              asymmetric_colorings(dihedral_grid(n)) // 8,
              cycle_index(dihedral_grid(n)))
