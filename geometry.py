"""Geometry for the arc-drawing puzzle.

Coordinate system: (row, col) with row 0 at the top.

Cell corners and edges:

    NW  --N--  NE
     |        |
     W        E
     |        |
    SW  --S--  SE

An arc is a quarter-circle of radius 1 centered at one of the 4 corners.
Its endpoints are the two corners adjacent to the center.

A cell with an arc is split into two sub-cells:
  - 'near' side: the quarter-disk (area pi/4), containing the center corner.
  - 'far'  side: the complement (area 1 - pi/4), containing the corner opposite the center.

A cell with no arc has a single sub-cell of area 1.
"""

from math import pi

# ---- Arc state encoding ----------------------------------------------------

# We use integers 0..4 for arc state. 0 = no arc; 1..4 = arc centered at the
# corresponding corner.
NO_ARC = 0
ARC_NW = 1
ARC_NE = 2
ARC_SW = 3
ARC_SE = 4

ARC_STATES = (NO_ARC, ARC_NW, ARC_NE, ARC_SW, ARC_SE)
ARCED_STATES = (ARC_NW, ARC_NE, ARC_SW, ARC_SE)

# Edges of a cell.
N, E, S, W = 'N', 'E', 'S', 'W'
EDGES = (N, E, S, W)

# Corners of a cell.
NW, NE, SW, SE = 'NW', 'NE', 'SW', 'SE'

# Map ARC_* state -> the corner where the arc is centered.
ARC_CENTER = {
    ARC_NW: NW,
    ARC_NE: NE,
    ARC_SW: SW,
    ARC_SE: SE,
}

# For each corner C, the two cell edges adjacent to C (i.e. the edges that
# meet at C). These are the edges bounding the 'near' side when the arc is
# centered at C.
EDGES_AT_CORNER = {
    NW: frozenset({N, W}),
    NE: frozenset({N, E}),
    SW: frozenset({S, W}),
    SE: frozenset({S, E}),
}

ALL_EDGES = frozenset(EDGES)

# ---- Sub-cell description --------------------------------------------------

# A "sub-cell" is identified by (row, col, side) where side is:
#   None   -- the whole cell (when there is no arc)
#   'near' -- the quarter-disk side (area pi/4)
#   'far'  -- the complement side  (area 1 - pi/4)
NEAR = 'near'
FAR = 'far'


def subcells(state):
    """Return the list of sides for a cell in the given arc state."""
    if state == NO_ARC:
        return [None]
    return [NEAR, FAR]


def subcell_edges(state, side):
    """Set of cell edges (subset of {N,E,S,W}) that the given sub-cell touches.

    Used for adjacency: two sub-cells in neighboring cells are in the same
    region iff both touch the shared edge.
    """
    if state == NO_ARC:
        assert side is None
        return ALL_EDGES
    center = ARC_CENTER[state]
    near_edges = EDGES_AT_CORNER[center]
    if side == NEAR:
        return near_edges
    if side == FAR:
        return ALL_EDGES - near_edges
    raise ValueError(f'bad side {side!r} for state {state}')


def subcell_area(state, side):
    if state == NO_ARC:
        return 1.0
    return (pi / 4) if side == NEAR else (1.0 - pi / 4)


def subcell_pi_count(state, side):
    """Coefficient of pi/4 in the sub-cell's area.

    Each cell contributes its area as (rational) + (pi_count) * pi/4.
    For a region's total area to be an integer, the sum of pi_counts must
    be 0, AND the rational part must be an integer.
    """
    if state == NO_ARC:
        return 0
    return 1 if side == NEAR else -1  # near = pi/4; far = 1 - pi/4


def subcell_rational(state, side):
    """Rational part of the sub-cell's area."""
    if state == NO_ARC:
        return 1
    return 0 if side == NEAR else 1


# ---- Cell-edge orientation (for adjacency) --------------------------------

# Opposite edges (used to translate "my N edge" -> "neighbor's S edge").
OPPOSITE_EDGE = {N: S, S: N, E: W, W: E}

# Direction deltas: moving across edge X takes you to the cell at offset...
EDGE_DELTA = {N: (-1, 0), S: (+1, 0), E: (0, +1), W: (0, -1)}
