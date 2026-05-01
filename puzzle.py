"""Puzzle data for the arc-drawing puzzle.

Grid is N x N. Coordinates are (row, col), row 0 = top, col 0 = left.

GREEN cells cannot contain arcs.
CLUES[(r, c)] = label  means the region containing >=half of cell (r,c) has score = label.

NOTE: Read directly from the puzzle image. If anything looks off, edit here.
"""

N = 9

# Green cells (no arc allowed). Clue cells can also be green.
GREEN = frozenset({
    (0, 0), (0, 3), (0, 5), (0, 7),
    (1, 3), (1, 6), (1, 8),
    (2, 0), (2, 6),
    (4, 4),
    (5, 7), (5, 8),
    (7, 0),
    (8, 1), (8, 4), (8, 8),
    # The "25" at (4,0) and "35" at (8,5) appear in green text in the image.
    (4, 0),
    (8, 5),
})

# Clue cells: (row, col) -> score of the region containing >=half of this cell.
CLUES = {
    (0, 2): 21,
    (1, 0): 21,
    (1, 4): 27,
    (1, 7): 25,
    (2, 1): 27,
    (2, 5): 15,
    (2, 8): 9,
    (4, 0): 25,
    (4, 3): 27,
    (4, 5): 45,
    (4, 8): 9,
    (6, 0): 9,
    (6, 3): 63,
    (6, 7): 45,
    (7, 1): 63,
    (7, 4): 9,
    (7, 8): 288,
    (8, 5): 35,
}


# --- Tiny example puzzle from the prompt (4x4), used as a sanity test. ---
# example grid:
#   3 . 9 .
#   . . . 6
#   8 . . .
#   . 6 . 24
# Greens (highlighted): (0,1), (0,3), (2,1), (2,3)
EXAMPLE_N = 4
EXAMPLE_GREEN = frozenset({(0, 1), (0, 3), (2, 1), (2, 3)})
EXAMPLE_CLUES = {
    (0, 0): 3,
    (0, 2): 9,
    (1, 3): 6,
    (2, 0): 8,
    (3, 1): 6,
    (3, 3): 24,
}
EXAMPLE_ANSWER = 18928
