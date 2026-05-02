"""Score computation: per-region score, per-cell score, final answer."""

from math import pi
from regions import find_regions, region_area
from boundary import collect_region_segments, count_smooth_pieces
from geometry import NO_ARC, NEAR, FAR, subcell_area


def region_score(states, region_subcells, n):
    """(score, integer_area, smooth_pieces) for a region.
    score = integer_area * smooth_pieces (or None if area is not integer).
    """
    _, is_int, int_area = region_area(states, region_subcells)
    if not is_int:
        return None, None, None
    segs = collect_region_segments(states, region_subcells, n)
    sp = count_smooth_pieces(segs)
    return int_area * sp, int_area, sp


def cell_main_subcell(state, side_choices):
    """Of the sub-cells of a cell, return the one that contains >= half the
    cell area. With our arc geometry, the 'near' side has area pi/4 ~ 0.785
    so it always wins over the 'far' side (area 1 - pi/4 ~ 0.215). For an
    unarced cell, the only sub-cell is the whole cell.
    """
    if state == NO_ARC:
        assert side_choices == [None]
        return None
    return NEAR


def evaluate(states, clues=None, *, strict=True):
    """Run the full scoring pipeline.

    Returns dict with:
      'regions'     : list of (region_id, [(r,c,side), ...])
      'region_score': dict region_id -> (score, area, pieces)
      'cell_score'  : NxN array of int (cell -> score of its >=half region)
      'row_sums'    : list of length n
      'col_sums'    : list of length n
      'answer'      : sum_of_squares(row_sums) + sum_of_squares(col_sums)
      'valid'       : bool, True iff all regions have integer area, no danglers,
                      and all clues are satisfied.
      'errors'      : list of strings explaining any constraint violations.

    If strict and !valid, raises with the errors.
    """
    n = len(states)
    sub2reg, regions = find_regions(states)
    errors = []

    # No-dangling: arced cells must have their two sides in different regions.
    for r in range(n):
        for c in range(n):
            if states[r][c] == NO_ARC:
                continue
            r_near = sub2reg[(r, c, NEAR)]
            r_far = sub2reg[(r, c, FAR)]
            if r_near == r_far:
                errors.append(f'dangling arc at ({r},{c}): both sides in region {r_near}')

    # Per-region scores; integer-area constraint.
    rscore = {}
    for rid, sc in regions.items():
        score, area, pieces = region_score(states, sc, n)
        rscore[rid] = (score, area, pieces)
        if score is None:
            errors.append(f'region {rid} has non-integer area (pi/4 not balanced)')

    # Cell -> region (the >=half side) -> score.
    cell_score = [[None] * n for _ in range(n)]
    for r in range(n):
        for c in range(n):
            st = states[r][c]
            sides = [None] if st == NO_ARC else [NEAR, FAR]
            main = cell_main_subcell(st, sides)
            rid = sub2reg[(r, c, main)]
            cell_score[r][c] = rscore[rid][0]

    # Clues.
    if clues:
        for (r, c), label in clues.items():
            got = cell_score[r][c]
            if got != label:
                errors.append(f'clue ({r},{c})={label} but got {got}')

    if strict and errors:
        raise ValueError('invalid configuration: ' + '; '.join(errors))

    row_sums = [sum(row) if all(x is not None for x in row) else None for row in cell_score]
    col_sums = [
        (sum(cell_score[r][c] for r in range(n))
         if all(cell_score[r][c] is not None for r in range(n))
         else None)
        for c in range(n)
    ]
    if any(x is None for x in row_sums) or any(x is None for x in col_sums):
        answer = None
    else:
        answer = sum(s * s for s in row_sums) + sum(s * s for s in col_sums)

    return {
        'regions': regions,
        'region_score': rscore,
        'cell_score': cell_score,
        'row_sums': row_sums,
        'col_sums': col_sums,
        'answer': answer,
        'valid': not errors,
        'errors': errors,
    }
