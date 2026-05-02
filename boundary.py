"""Boundary tracing and smooth-piece counting for regions.

Coordinates: cell (r, c) occupies the unit square with corners
    NW = (r,   c)
    NE = (r,   c+1)
    SW = (r+1, c)
    SE = (r+1, c+1)
in (row, col) space (row increasing downward).

A region's boundary is composed of:
  * arc segments (each arc separates two regions),
  * outer-grid edge segments (cell edges that border the outside of the grid
    AND lie on the side of the region's sub-cell).

Cell-internal edges between adjacent cells are NEVER region boundaries -- the
two sub-cells across them are always in the same region (we only get boundary
crossings at arcs, and arcs only touch cell corners).

Number of smooth pieces in a region:
  Walk each closed boundary loop. Count vertices where the tangent direction
  changes (corners). For a loop with k corners (k >= 1), the loop contributes
  k smooth pieces. For a corner-free loop (purely smooth, e.g. a circle made
  of 4 arcs), the loop contributes 1 smooth piece.
  Total smooth pieces = sum over all loops.
"""

from geometry import (
    NO_ARC, NEAR, FAR,
    N, E, S, W,
    NW, NE, SW, SE,
    ARC_CENTER, EDGES_AT_CORNER,
    subcells, subcell_edges,
)


# --- Geometry helpers -------------------------------------------------------

# Cell corners as (row, col) lattice points. For cell (r,c):
def corner_point(r, c, corner):
    if corner == NW: return (r, c)
    if corner == NE: return (r, c + 1)
    if corner == SW: return (r + 1, c)
    if corner == SE: return (r + 1, c + 1)
    raise ValueError(corner)


# Endpoints of the two corners of an edge.
EDGE_ENDPOINTS = {
    N: (NW, NE),
    S: (SW, SE),
    E: (NE, SE),
    W: (NW, SW),
}


# For an arc centered at corner C, its two endpoint corners.
def arc_endpoints(center):
    if center == NW: return (NE, SW)
    if center == NE: return (NW, SE)
    if center == SW: return (NW, SE)
    if center == SE: return (NE, SW)
    raise ValueError(center)


# Tangent direction (unit-ish vector in (drow, dcol)) of an arc at one of its
# endpoints, oriented "going INTO the arc, away from the endpoint".
#
# For arc centered at corner C with endpoint P (an adjacent corner of the cell):
# the tangent at P is perpendicular to (P - C). It points "along the cell edge
# from P that goes away from C".
#
# Concretely:
#   center NW, endpoint NE: P-C points east -> tangent is along the E edge
#     leaving NE southward, i.e. (drow, dcol) = (+1, 0).
#   center NW, endpoint SW: P-C points south -> tangent leaves SW eastward, (0, +1).
#   center NE, endpoint NW: P-C points west -> tangent leaves NW southward, (+1, 0).
#   center NE, endpoint SE: P-C points south -> tangent leaves SE westward, (0, -1).
#   center SW, endpoint NW: P-C points north -> tangent leaves NW eastward, (0, +1).
#   center SW, endpoint SE: P-C points east -> tangent leaves SE northward, (-1, 0).
#   center SE, endpoint NE: P-C points north -> tangent leaves NE westward, (0, -1).
#   center SE, endpoint SW: P-C points west -> tangent leaves SW northward, (-1, 0).
ARC_TANGENT_INTO_ARC = {
    (NW, NE): (+1,  0),
    (NW, SW): ( 0, +1),
    (NE, NW): (+1,  0),
    (NE, SE): ( 0, -1),
    (SW, NW): ( 0, +1),
    (SW, SE): (-1,  0),
    (SE, NE): ( 0, -1),
    (SE, SW): (-1,  0),
}


# --- Boundary segments ------------------------------------------------------

# Each region-boundary segment is one of:
#   ('arc', r, c, center)      -- the arc in cell (r,c) centered at `center`
#   ('edge', r, c, edge)       -- a cell edge of cell (r,c) on the outer grid boundary
#
# Each arc segment is shared between TWO sub-cells (the near and far sides of
# the arced cell), but those two sub-cells are in DIFFERENT regions (assuming
# no dangling). So an arc appears on the boundary of two regions.
#
# Each outer grid-edge segment belongs to exactly ONE region (the region of
# the sub-cell of cell (r,c) that touches `edge`).


def collect_region_segments(states, region_subcells, n):
    """For a single region (list of (r, c, side) sub-cells), return the list of
    boundary segments.

    Each segment is a tuple as described above, plus we record both endpoints
    (corner labels) and the tangent vectors at each endpoint so that the
    smooth-piece counter doesn't have to recompute them.
    """
    segs = []
    in_region = set(region_subcells)

    # Collect arcs first: for each cell, if it has an arc and either side is
    # in this region, the arc is a boundary segment of this region.
    seen_arcs = set()
    for (r, c, side) in region_subcells:
        st = states[r][c]
        if st == NO_ARC:
            continue
        if (r, c) in seen_arcs:
            continue
        seen_arcs.add((r, c))
        center = ARC_CENTER[st]
        ep1, ep2 = arc_endpoints(center)
        # Tangents at the two endpoints, both pointing INTO the arc.
        t1 = ARC_TANGENT_INTO_ARC[(center, ep1)]
        t2 = ARC_TANGENT_INTO_ARC[(center, ep2)]
        p1 = corner_point(r, c, ep1)
        p2 = corner_point(r, c, ep2)
        segs.append({
            'kind': 'arc',
            'cell': (r, c),
            'center': center,
            'ep1': p1, 't1': t1,   # tangent at p1 going INTO arc (toward p2)
            'ep2': p2, 't2': t2,   # tangent at p2 going INTO arc (toward p1)
        })

    # Outer-grid edges that border this region.
    for (r, c, side) in region_subcells:
        st = states[r][c]
        # Which edges does this sub-cell touch?
        edges_touched = subcell_edges(st, side)
        for edge in edges_touched:
            # Is `edge` on the outer boundary of the grid?
            if edge == N and r == 0:
                pass
            elif edge == S and r == n - 1:
                pass
            elif edge == W and c == 0:
                pass
            elif edge == E and c == n - 1:
                pass
            else:
                continue
            # Add this outer-edge segment.
            ep_a, ep_b = EDGE_ENDPOINTS[edge]
            pa = corner_point(r, c, ep_a)
            pb = corner_point(r, c, ep_b)
            # Tangent: edges are straight, tangent is constant along them.
            # For corner-counting we only care about direction along the edge.
            # Make tangent point from pa to pb.
            tangent = (pb[0] - pa[0], pb[1] - pa[1])
            # Normalize to unit step:
            mag = abs(tangent[0]) + abs(tangent[1])
            tangent = (tangent[0] // mag, tangent[1] // mag)
            segs.append({
                'kind': 'edge',
                'cell': (r, c),
                'edge': edge,
                'ep1': pa, 't1': tangent,         # tangent leaving pa toward pb
                'ep2': pb, 't2': (-tangent[0], -tangent[1]),  # leaving pb toward pa
            })
    return segs


def count_smooth_pieces(segments):
    """Given the list of boundary segments for ONE region, return the number
    of smooth pieces.

    Algorithm: build an undirected multigraph on corner-points. Each segment
    contributes one edge between ep1 and ep2. At each vertex, pair up the
    incident segment-ends so that paired ends have OPPOSITE tangent vectors
    (so a smooth walk passes through the pair). Then each connected pair-chain
    is one smooth piece (or a closed smooth loop counts as 1 piece).

    We assume the region's boundary is well-formed: at every vertex, the
    segment-ends can be paired up such that tangents face opposite directions
    (i.e. the boundary walk is consistent).
    """
    # Build per-vertex list of (segment_index, end_label) where end_label is 1 or 2.
    by_vertex = {}
    for i, seg in enumerate(segments):
        by_vertex.setdefault(seg['ep1'], []).append((i, 1))
        by_vertex.setdefault(seg['ep2'], []).append((i, 2))

    def end_tangent(seg, label):
        """Tangent at the end, pointing AWAY from the segment (i.e. outgoing).
        Stored t1 / t2 already point INTO the arc / along the edge from that
        endpoint, so they are the OUTGOING direction from the vertex."""
        return seg['t1'] if label == 1 else seg['t2']

    # At each vertex, pair ends with opposite outgoing tangents (so they form
    # a smooth pass-through). If a pair's tangents are NOT opposite, it's a
    # corner -- in that case the two ends still get paired (they're consecutive
    # around the vertex) but they belong to DIFFERENT smooth pieces.
    #
    # Approach: union-find over (segment_index, end_label) pairs. We union
    # ends at the same vertex that have opposite tangents. At the end, each
    # connected component of THE INCIDENCE-PAIR graph is one smooth piece.
    #
    # ...but we still need to follow ends through segments: each segment
    # connects its end1 to its end2 along itself (which is by definition
    # smooth -- arcs and edges are individually smooth). So union those too.

    parent = {}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for i in range(len(segments)):
        for end in (1, 2):
            parent[(i, end)] = (i, end)

    # Union the two ends of each segment (smooth along the segment itself).
    for i in range(len(segments)):
        union((i, 1), (i, 2))

    # At each vertex, pair ends with opposite outgoing tangents.
    for vertex, ends in by_vertex.items():
        # Group by tangent direction.
        for j in range(len(ends)):
            ij, lj = ends[j]
            tj = end_tangent(segments[ij], lj)
            for k in range(j + 1, len(ends)):
                ik, lk = ends[k]
                tk = end_tangent(segments[ik], lk)
                if (tj[0] + tk[0], tj[1] + tk[1]) == (0, 0):
                    # Opposite tangents -> smooth pass-through.
                    union((ij, lj), (ik, lk))

    # Count distinct components.
    roots = set()
    for i in range(len(segments)):
        roots.add(find((i, 1)))
    return len(roots)
