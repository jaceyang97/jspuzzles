"""Backtracking solver for the arc-drawing puzzle.

Performance ideas in this revision:
  * Skip smooth-piece counting for regions with no clue cell. Only the
    integer-area check (very cheap) is done unless a clue is present.
  * Keep an incremental Union-Find with rollback so each cell assignment
    only does O(1) work plus a few unions.
  * Track "open edges" per region (count of edges that leak to unassigned
    cells). When this hits zero, the region is closed; validate immediately.
"""

import sys
from collections import defaultdict
from geometry import (
    NO_ARC, NEAR, FAR, ARC_STATES, N, E, S, W,
    EDGE_DELTA, OPPOSITE_EDGE,
    subcells, subcell_edges, subcell_pi_count, subcell_rational,
)
from boundary import collect_region_segments, count_smooth_pieces


# ---- DSU with rollback -----------------------------------------------------

class RollbackDSU:
    """Union-find with explicit log so unions can be undone."""

    def __init__(self):
        self.parent = {}
        self.rank = {}
        self.log = []   # list of operations to undo

    def add(self, x):
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0
            self.log.append(('add', x))

    def find(self, x):
        # No path compression (would complicate rollback).
        while self.parent[x] != x:
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            self.log.append(('noop',))
            return False
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        old_parent_rb = self.parent[rb]
        old_rank_ra = self.rank[ra]
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1
        self.log.append(('union', rb, old_parent_rb, ra, old_rank_ra))
        return True

    def checkpoint(self):
        return len(self.log)

    def rollback(self, cp):
        while len(self.log) > cp:
            op = self.log.pop()
            if op[0] == 'add':
                _, x = op
                del self.parent[x]
                del self.rank[x]
            elif op[0] == 'noop':
                pass
            elif op[0] == 'union':
                _, rb, old_parent_rb, ra, old_rank_ra = op
                self.parent[rb] = old_parent_rb
                self.rank[ra] = old_rank_ra


# ---- Solver ----------------------------------------------------------------

class Solver:
    def __init__(self, n, green_cells, clues, *, verbose=False):
        self.n = n
        self.green = set(green_cells)
        self.clues = dict(clues)
        self.verbose = verbose

    def domain(self, r, c):
        if (r, c) in self.green:
            return [NO_ARC]
        return list(ARC_STATES)

    def solve(self, max_solutions=1, max_steps=None):
        n = self.n
        states = [[NO_ARC] * n for _ in range(n)]
        order = [(r, c) for r in range(n) for c in range(n)]

        # Per-region stats keyed by DSU root: (rational_sum, pi_quarter_sum,
        # near_cells_set, far_cells_set, open_edges_count, clue_cells_in_region)
        # We track them via a single structure indexed by root, recomputed by
        # re-scan each step (still cheap relative to find_regions full).

        self.steps = 0
        self.solutions = []
        self._max_steps = max_steps

        # We'll keep a DSU over assigned sub-cells.
        self.dsu = RollbackDSU()

        try:
            self._backtrack(states, order, 0, max_solutions)
        except _StopSearch:
            pass
        return self.solutions

    # Build a list of (subcell, edges_touched) for a cell+state.
    def _cell_subcells(self, st):
        out = []
        for side in subcells(st):
            out.append((side, subcell_edges(st, side)))
        return out

    def _backtrack(self, states, order, idx, max_solutions):
        if len(self.solutions) >= max_solutions:
            return
        self.steps += 1
        if self._max_steps is not None and self.steps > self._max_steps:
            raise _StopSearch()
        if self.verbose and self.steps % 100000 == 0:
            print(f'  steps={self.steps:,} idx={idx}/{len(order)}', file=sys.stderr, flush=True)

        if idx == len(order):
            ok, info = self._validate_full(states)
            if ok:
                self.solutions.append(([row[:] for row in states], info))
            return

        r, c = order[idx]
        n = self.n
        domain = self.domain(r, c)
        for st in domain:
            cp = self.dsu.checkpoint()
            states[r][c] = st
            ok = self._add_cell_and_check(states, r, c, idx, order)
            if ok:
                self._backtrack(states, order, idx + 1, max_solutions)
                if len(self.solutions) >= max_solutions:
                    return
            states[r][c] = NO_ARC
            self.dsu.rollback(cp)

    def _add_cell_and_check(self, states, r, c, idx, order):
        """Add cell (r,c) to the DSU, do unions, then validate any newly-closed
        regions. Returns True iff still consistent.
        """
        n = self.n
        st = states[r][c]
        # Add sub-cells.
        for side in subcells(st):
            self.dsu.add((r, c, side))
        # Unions with already-assigned neighbors via N and W (cells assigned earlier in row-major).
        # Also need to union via S and E if those neighbors are assigned (but in
        # row-major they're not yet). We're row-major; only N and W are
        # already-assigned at this point (possibly).
        for edge in (N, W, S, E):
            dr, dc = EDGE_DELTA[edge]
            nr, nc = r + dr, c + dc
            if not (0 <= nr < n and 0 <= nc < n):
                continue
            # Is the neighbor cell already assigned? In row-major order, cells
            # before (r, c) in the order list are assigned. We need to know
            # which cells are assigned. Use the index argument.
            if not self._is_assigned(nr, nc, idx, order):
                continue
            nst = states[nr][nc]
            nedge = OPPOSITE_EDGE[edge]
            for side in subcells(st):
                if edge in subcell_edges(st, side):
                    for nside in subcells(nst):
                        if nedge in subcell_edges(nst, nside):
                            self.dsu.union((r, c, side), (nr, nc, nside))

        # Validate any newly-closed regions.
        # Build per-root stats by scanning all assigned subcells. This is O(n^2)
        # per step but fast in pure Python relative to find_regions full.
        return self._validate_partial(states, idx, order)

    def _is_assigned(self, r, c, idx, order):
        # Cells with index < idx in `order` are assigned. We need quick lookup;
        # build a position map lazily.
        pos = self.__dict__.setdefault('_pos', {})
        if not pos:
            for i, (rr, cc) in enumerate(order):
                pos[(rr, cc)] = i
        return pos[(r, c)] < idx

    def _validate_partial(self, states, idx, order):
        n = self.n
        # Determine open-edge count per root: count of edges of assigned cells
        # that border UNASSIGNED neighbors and whose subcell touches that edge.
        # Per-root accumulators:
        rat = defaultdict(int)
        pi_q = defaultdict(int)
        cells_with_near = defaultdict(set)
        cells_with_far = defaultdict(set)
        open_edges = defaultdict(int)
        members = defaultdict(list)
        clues_in = defaultdict(list)

        # Scan assigned cells.
        for i in range(idx + 1):
            r, c = order[i]
            st = states[r][c]
            for side in subcells(st):
                root = self.dsu.find((r, c, side))
                rat[root] += subcell_rational(st, side)
                pi_q[root] += subcell_pi_count(st, side)
                if side == NEAR:
                    cells_with_near[root].add((r, c))
                elif side == FAR:
                    cells_with_far[root].add((r, c))
                members[root].append((r, c, side))
                # Open edges: this subcell's edges that go to unassigned cells.
                edges_touched = subcell_edges(st, side)
                for edge in edges_touched:
                    dr, dc = EDGE_DELTA[edge]
                    nr, nc = r + dr, c + dc
                    if not (0 <= nr < n and 0 <= nc < n):
                        continue  # outer boundary, not a leak
                    if not self._is_assigned(nr, nc, idx + 1, order):
                        open_edges[root] += 1
            if (r, c) in self.clues:
                main_side = None if st == NO_ARC else NEAR
                root = self.dsu.find((r, c, main_side))
                clues_in[root].append((r, c, self.clues[(r, c)]))

        # Open-region prunes (cheap, apply to all regions):
        for root in list(rat.keys()):
            # No-dangling: if any arced cell has BOTH sides in this region,
            # it's stuck (DSU only grows), so fail now.
            if cells_with_near[root] & cells_with_far[root]:
                return False
            # Clue area bound: clue cell's region's eventual area >= current
            # rat (since adding cells only adds non-negative rat). So if
            # rat > label, no way: pieces * area = label, area >= rat,
            # so pieces <= label/rat < 1 -> impossible (pieces >= 1).
            for r, c, label in clues_in[root]:
                if rat[root] > label:
                    return False

        # Validate closed regions.
        for root in list(rat.keys()):
            if open_edges[root] > 0:
                continue
            # Closed region.
            # Integer area.
            if pi_q[root] != 0:
                return False
            area = rat[root]
            if clues_in[root]:
                segs = collect_region_segments(states, members[root], n)
                pieces = count_smooth_pieces(segs)
                for r, c, label in clues_in[root]:
                    if area * pieces != label:
                        return False

        return True

    def _validate_full(self, states):
        from regions import find_regions
        from score import region_score
        n = self.n
        sub2reg, regions = find_regions(states)
        for r in range(n):
            for c in range(n):
                st = states[r][c]
                if st == NO_ARC: continue
                if sub2reg[(r, c, NEAR)] == sub2reg[(r, c, FAR)]:
                    return False, None
        rscore = {}
        for rid, sc in regions.items():
            score, area, pieces = region_score(states, sc, n)
            if score is None:
                return False, None
            rscore[rid] = (score, area, pieces)
        for (r, c), label in self.clues.items():
            st = states[r][c]
            main = None if st == NO_ARC else NEAR
            if rscore[sub2reg[(r, c, main)]][0] != label:
                return False, None
        return True, {'regions': regions, 'sub2reg': sub2reg, 'region_score': rscore}


class _StopSearch(Exception):
    pass
