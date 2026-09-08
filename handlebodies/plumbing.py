"""Trees whose plumbing of cylinders realises a given surface.

A surface is built by plumbing copies of DT*S^1 -- cylinders -- one per vertex
of a graph, glued where its edges say.  Such a plumbing has Euler
characteristic minus the number of edges, and for a tree the cores are a basis
of H_1 whose intersection form is the tree's skew adjacency matrix.  Genus is
half the rank of that form and the ends are its corank plus one, and the rank
of a forest's skew adjacency matrix is twice its maximum matching.  So

    genus = maximum matching,      ends = vertices - 2 * matching + 1,

which inverts to: the trees realising a surface of genus g with n ends are
those on 2g + n - 1 vertices whose maximum matching is g.  Equivalently, the
matching must miss exactly n - 1 vertices.

The surface does not determine the tree -- a skeleton is a choice, not an
invariant -- so there is generally a family of them, and in dimension two they
all present the same Liouville surface.  Two cases are rigid: a tree with no
two disjoint edges is a star, so genus one is always the star on n + 1
vertices; and no tree at all has genus zero except the single cylinder.
"""

import itertools
from functools import lru_cache

#: Beyond this the trees are more than are worth listing, and the page asks
#: for them on every computation whether or not anyone opens the section.
#: Twelve vertices is 551 trees and a tenth of a second; fourteen is seven
#: times that.
VERTEX_LIMIT = 12


def _nodes(tree):
    return 1 + sum(_nodes(child) for child in tree)


@lru_cache(maxsize=None)
def _rooted(n):
    """Rooted trees on `n` nodes, each a tuple of its children's trees."""
    if n == 1:
        return ((),)
    return tuple(_forests(n - 1, n - 1))


@lru_cache(maxsize=None)
def _forests(total, largest):
    """Multisets of rooted trees on `total` nodes, none larger than `largest`.

    Kept in nonincreasing order so that each multiset is produced once rather
    than once per permutation of its parts.
    """
    if total == 0:
        return ((),)
    out = []
    for size in range(min(total, largest), 0, -1):
        for tree in _rooted(size):
            for rest in _forests(total - size, size):
                if rest and (_nodes(rest[0]), rest[0]) > (size, tree):
                    continue
                out.append((tree,) + rest)
    return tuple(out)


def _edges(tree):
    """A rooted tree as an edge list on 0..n-1, the root being 0."""
    edges, nxt = [], 1

    def walk(node, index):
        nonlocal nxt
        for child in node:
            here = nxt
            nxt += 1
            edges.append((index, here))
            walk(child, here)

    walk(tree, 0)
    return edges


def _canonical(edges, n):
    """A free tree's canonical form: rooted at its centre, children sorted."""
    adjacent = [[] for _ in range(n)]
    for a, b in edges:
        adjacent[a].append(b)
        adjacent[b].append(a)
    degree = [len(a) for a in adjacent]
    gone, left = set(), n
    layer = [v for v in range(n) if degree[v] <= 1]
    while left > 2:
        nxt = []
        for v in layer:
            gone.add(v)
            left -= 1
            for w in adjacent[v]:
                if w not in gone:
                    degree[w] -= 1
                    if degree[w] == 1:
                        nxt.append(w)
        layer = nxt
    centres = [v for v in range(n) if v not in gone]

    def code(v, parent):
        return '(' + ''.join(sorted(code(w, v) for w in adjacent[v]
                                    if w != parent)) + ')'

    if len(centres) == 1:
        return code(centres[0], -1)
    a, b = centres
    return min(code(a, b) + code(b, a), code(b, a) + code(a, b))


def matching_number(edges, n):
    """The largest set of disjoint edges, which for a tree is a greedy count.

    Matching leaves to their parents first is optimal on a tree: a leaf's only
    edge is the one to its parent, so nothing is lost by taking it.
    """
    if n < 2:
        return 0
    adjacent = [[] for _ in range(n)]
    for a, b in edges:
        adjacent[a].append(b)
        adjacent[b].append(a)
    parent, order, stack, seen = [-1] * n, [], [0], [False] * n
    seen[0] = True
    while stack:
        v = stack.pop()
        order.append(v)
        for w in adjacent[v]:
            if not seen[w]:
                seen[w] = True
                parent[w] = v
                stack.append(w)
    used, size = [False] * n, 0
    for v in reversed(order):                     # deepest first
        p = parent[v]
        if p >= 0 and not used[v] and not used[p]:
            used[v] = used[p] = True
            size += 1
    return size


def free_trees(n):
    """Every tree on `n` vertices, once each, as an edge list."""
    if n <= 0:
        return []
    seen = {}
    for tree in _rooted(n):
        edges = _edges(tree)
        seen.setdefault(_canonical(edges, n), edges)
    return list(seen.values())


def _compositions(total, parts):
    """Ways of splitting `total` into `parts` positive pieces."""
    for cuts in itertools.combinations(range(1, total), parts - 1):
        edges = (0,) + cuts + (total,)
        yield tuple(edges[i + 1] - edges[i] for i in range(parts))


def some_realising(genus, ends, wanted):
    """A few trees realising the surface, built rather than searched for.

    Enumerating every tree on many vertices is the expensive part, so these are
    constructed instead: take any tree of `genus` centres and hang the
    remaining vertices off them as leaves, every centre carrying at least one.
    The centres are then a vertex cover, so the matching is at most their
    number, and a leaf apiece gives that many disjoint edges, so it is exactly
    their number.

    Not every tree has this shape -- A5 does not, its cover being two vertices
    that are not adjacent and its middle vertex being no leaf -- so this is a
    supply of examples, not a census.
    """
    if genus < 1:
        return []
    leaves = genus + ends - 1
    if leaves < genus:
        return []
    found = {}
    for centres in free_trees(genus):
        for share in _compositions(leaves, genus):
            edges, nxt = list(centres), genus
            for centre, count in enumerate(share):
                for _ in range(count):
                    edges.append((centre, nxt))
                    nxt += 1
            found.setdefault(_canonical(edges, nxt), edges)
            if len(found) >= wanted:
                return list(found.values())
    return list(found.values())


def realising(genus, ends, sample=20):
    """The trees whose plumbing gives a surface of this genus and this many ends.

    Complete while the vertices are few enough to enumerate them all, and past
    that a handful built by hand: knowing one tree that works is worth more
    than knowing that there are too many to count.

    Genus one is complete either way.  A tree with no two disjoint edges is a
    star, so there is only ever the one.
    """
    n = 2 * genus + ends - 1
    if n < 1:
        return {'vertices': n, 'trees': [], 'listed': True}
    if n <= VERTEX_LIMIT:
        trees = [e for e in free_trees(n) if matching_number(e, n) == genus]
        return {'vertices': n, 'trees': trees, 'listed': True}
    trees = some_realising(genus, ends, sample)
    return {'vertices': n, 'trees': trees, 'listed': genus == 1}
