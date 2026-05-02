"""Region detection for the arc-drawing puzzle.

A region is a maximal set of sub-cells reachable from one another by crossing
shared cell edges. Arcs only touch at cell corners, so a shared cell edge is
fully passable -- two sub-cells across it are in the same region iff each
touches that shared edge.
"""

from geometry import (
    NO_ARC, NEAR, FAR,
    N, E, S, W,
    EDGES, OPPOSITE_EDGE, EDGE_DELTA,
    subcells, subcell_edges, subcell_pi_count, subcell_rational,
)


class DSU:
    def __init__(self):
        self.parent = {}

    def add(self, x):
        if x not in self.parent:
            self.parent[x] = x

    def find(self, x):
        # Path compression.
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[x] != root:
            self.parent[x], x = root, self.parent[x]
        return root

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


def find_regions(states):
    """Given an NxN list-of-lists of arc states, return:

      sub2region : dict (r, c, side) -> region_id (small integer)
      regions    : dict region_id -> list of (r, c, side) sub-cells

    Where `side` is None for unarced cells, or 'near'/'far' for arced cells.
    """
    n = len(states)
    dsu = DSU()

    # Add every sub-cell.
    for r in range(n):
        for c in range(n):
            for side in subcells(states[r][c]):
                dsu.add((r, c, side))

    # Union sub-cells across shared cell edges.
    # Iterate each cell; for each of its edges (only need to do N and W to
    # avoid double-counting), connect to the neighbor.
    for r in range(n):
        for c in range(n):
            st = states[r][c]
            for edge in (N, W):  # neighbors above and to the left (avoid double-count)
                dr, dc = EDGE_DELTA[edge]
                nr, nc = r + dr, c + dc
                if not (0 <= nr < n and 0 <= nc < n):
                    continue
                neighbor_st = states[nr][nc]
                neighbor_edge = OPPOSITE_EDGE[edge]
                # Find sub-cells of (r,c) touching `edge` and (nr,nc) touching neighbor_edge.
                for side in subcells(st):
                    if edge in subcell_edges(st, side):
                        for nside in subcells(neighbor_st):
                            if neighbor_edge in subcell_edges(neighbor_st, nside):
                                dsu.union((r, c, side), (nr, nc, nside))

    # Materialize regions, with stable ids in row-major scan order.
    sub2region = {}
    regions = {}
    next_id = 0
    for r in range(n):
        for c in range(n):
            for side in subcells(states[r][c]):
                root = dsu.find((r, c, side))
                if root not in sub2region:
                    sub2region[root] = next_id
                    regions[next_id] = []
                    next_id += 1
                rid = sub2region[root]
                regions[rid].append((r, c, side))

    # Re-key sub2region from sub-cell -> region id.
    out = {}
    for r in range(n):
        for c in range(n):
            for side in subcells(states[r][c]):
                out[(r, c, side)] = sub2region[dsu.find((r, c, side))]
    return out, regions


def region_area(states, region):
    """Return (numeric_area, is_integer, integer_area_or_None) for a region.

    Total area = rational_sum + (pi_quarter_count) * pi/4. Integer iff
    pi_quarter_count == 0 (since pi is irrational and rational parts are
    always integers in this puzzle).
    """
    from math import pi as _pi
    rat = 0
    pq = 0
    for r, c, side in region:
        st = states[r][c]
        rat += subcell_rational(st, side)
        pq += subcell_pi_count(st, side)
    is_int = (pq == 0)
    return rat + pq * (_pi / 4), is_int, (rat if is_int else None)
