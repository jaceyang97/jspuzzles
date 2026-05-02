"""Backtracking solver with incremental region tracking.

Per-root state is maintained incrementally; rollback restores state on
backtrack. Each cell assignment does only the work proportional to that
cell's sub-cells (constant), not a full O(n^2) re-scan.
"""

import sys
from geometry import (
    NO_ARC, NEAR, FAR, ARC_STATES, N, E, S, W,
    ARC_NW, ARC_NE, ARC_SW, ARC_SE,
    EDGE_DELTA, OPPOSITE_EDGE,
    subcells, subcell_edges, subcell_pi_count, subcell_rational,
)
from boundary import collect_region_segments, count_smooth_pieces


# Arc states whose endpoint lies at each corner of its cell.
_STATES_WITH_ENDPOINT_AT = {
    'NW': frozenset({ARC_NE, ARC_SW}),
    'NE': frozenset({ARC_NW, ARC_SE}),
    'SW': frozenset({ARC_NW, ARC_SE}),
    'SE': frozenset({ARC_NE, ARC_SW}),
}


def _endpoint_count_at_vertex(states, r, c, n):
    cnt = 0
    if r - 1 >= 0 and c - 1 >= 0 and states[r - 1][c - 1] in _STATES_WITH_ENDPOINT_AT['SE']: cnt += 1
    if r - 1 >= 0 and c < n     and states[r - 1][c]     in _STATES_WITH_ENDPOINT_AT['SW']: cnt += 1
    if r < n     and c - 1 >= 0 and states[r][c - 1]     in _STATES_WITH_ENDPOINT_AT['NE']: cnt += 1
    if r < n     and c < n      and states[r][c]         in _STATES_WITH_ENDPOINT_AT['NW']: cnt += 1
    return cnt


# ---- Incremental DSU + per-root stats with rollback log -------------------

class _Log:
    __slots__ = ('items',)
    def __init__(self):
        self.items = []
    def append(self, op):
        self.items.append(op)
    def checkpoint(self):
        return len(self.items)
    def restore(self, cp, dsu):
        items = self.items
        while len(items) > cp:
            op = items.pop()
            dsu._undo(op)


class IncrementalDSU:
    def __init__(self):
        self.parent = {}
        self.rank = {}
        # Per-root stats:
        self.rat = {}
        self.pi_q = {}
        self.near_cells = {}    # set
        self.far_cells = {}     # set
        self.clues = {}         # list of (r, c, label)
        self.open_edges = {}    # int
        self.members = {}       # list of (r, c, side) -- for boundary tracing
        self.log = _Log()

    # ---- low-level ops ----
    def _set(self, name, key, val):
        d = getattr(self, name)
        old = d.get(key, _MISSING)
        d[key] = val
        self.log.append(('set', name, key, old))

    def _del(self, name, key):
        d = getattr(self, name)
        old = d[key]
        del d[key]
        self.log.append(('del', name, key, old))

    def _undo(self, op):
        if op[0] == 'set':
            _, name, key, old = op
            d = getattr(self, name)
            if old is _MISSING:
                del d[key]
            else:
                d[key] = old
        elif op[0] == 'del':
            _, name, key, old = op
            d = getattr(self, name)
            d[key] = old

    def checkpoint(self):
        return self.log.checkpoint()

    def rollback(self, cp):
        self.log.restore(cp, self)

    # ---- DSU ----
    def find(self, x):
        # No path compression (would mutate parent without logging).
        while self.parent[x] != x:
            x = self.parent[x]
        return x

    def add(self, x):
        self._set('parent', x, x)
        self._set('rank', x, 0)

    def init_root_stats(self, x, *, rat, pi_q, is_near, is_far, cell, clue, init_open):
        """Initialize per-root stats for a brand-new singleton root x."""
        self._set('rat', x, rat)
        self._set('pi_q', x, pi_q)
        self._set('near_cells', x, {cell} if is_near else set())
        self._set('far_cells', x, {cell} if is_far else set())
        self._set('clues', x, list(clue) if clue else [])
        self._set('open_edges', x, init_open)
        self._set('members', x, [x])

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return ra
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        # Merge stats from rb -> ra.
        self._set('rat', ra, self.rat[ra] + self.rat[rb])
        self._set('pi_q', ra, self.pi_q[ra] + self.pi_q[rb])
        # Sets: shallow union. Log original ra's set so we can restore.
        new_near = self.near_cells[ra] | self.near_cells[rb]
        new_far = self.far_cells[ra] | self.far_cells[rb]
        self._set('near_cells', ra, new_near)
        self._set('far_cells', ra, new_far)
        self._set('clues', ra, self.clues[ra] + self.clues[rb])
        self._set('open_edges', ra, self.open_edges[ra] + self.open_edges[rb])
        self._set('members', ra, self.members[ra] + self.members[rb])
        # Drop rb's stats (free memory + ensures no orphans).
        self._del('rat', rb)
        self._del('pi_q', rb)
        self._del('near_cells', rb)
        self._del('far_cells', rb)
        self._del('clues', rb)
        self._del('open_edges', rb)
        self._del('members', rb)
        # Promote rb's parent.
        self._set('parent', rb, ra)
        if self.rank[ra] == self.rank[rb]:
            self._set('rank', ra, self.rank[ra] + 1)
        return ra

    def add_open_edges(self, x, delta):
        """Adjust open-edge count on the root of x by delta (can be negative)."""
        r = self.find(x)
        self._set('open_edges', r, self.open_edges[r] + delta)
        return r


_MISSING = object()


# ---- Solver ----------------------------------------------------------------

class Solver:
    def __init__(self, n, green_cells, clues, *, verbose=False, log_every=10000):
        self.n = n
        self.green = set(green_cells)
        self.clues = dict(clues)
        self.verbose = verbose
        self.log_every = log_every

    def domain(self, r, c):
        if (r, c) in self.green:
            return [NO_ARC]
        return list(ARC_STATES)

    def solve(self, max_solutions=1, max_steps=None):
        n = self.n
        states = [[NO_ARC] * n for _ in range(n)]
        order = [(r, c) for r in range(n) for c in range(n)]
        self.steps = 0
        self.solutions = []
        self._max_steps = max_steps
        self.dsu = IncrementalDSU()

        # Build a position map for "is assigned" queries.
        self._pos = {(r, c): i for i, (r, c) in enumerate(order)}

        try:
            self._backtrack(states, order, 0, max_solutions)
        except _StopSearch:
            pass
        return self.solutions

    def _is_assigned(self, r, c, idx):
        return self._pos[(r, c)] < idx

    def _backtrack(self, states, order, idx, max_solutions):
        if len(self.solutions) >= max_solutions:
            return
        self.steps += 1
        if self._max_steps is not None and self.steps > self._max_steps:
            raise _StopSearch()
        if self.verbose and self.steps % self.log_every == 0:
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
            ok = self._add_cell(states, r, c, idx)
            if ok:
                ok2 = self._after_add_check(states, r, c, idx)
                if ok2:
                    self._backtrack(states, order, idx + 1, max_solutions)
                    if len(self.solutions) >= max_solutions:
                        return
            states[r][c] = NO_ARC
            self.dsu.rollback(cp)

    # ---- Add a cell incrementally to the DSU + stats ----
    def _add_cell(self, states, r, c, idx):
        """Add cell (r, c) to DSU. Update open-edges of newly-closed neighbor
        edges. Return True if no immediate inconsistency.

        Inconsistency conditions raised here:
          * Newly-formed dangle: arced cell where DSU merge places NEAR and
            FAR in the same root.
        """
        n = self.n
        st = states[r][c]
        # Initialize each sub-cell as a singleton root.
        for side in subcells(st):
            self.dsu.add((r, c, side))
            edges_touched = subcell_edges(st, side)
            # Count open edges of this sub-cell (edges to unassigned neighbors).
            open_init = 0
            for edge in edges_touched:
                dr, dc = EDGE_DELTA[edge]
                nr, nc = r + dr, c + dc
                if not (0 <= nr < n and 0 <= nc < n):
                    continue  # outer boundary, not open
                if not self._is_assigned(nr, nc, idx):
                    open_init += 1
            # Initial stats for this singleton.
            clue = []
            main_side = None if st == NO_ARC else NEAR
            if (r, c) in self.clues and side == main_side:
                clue = [(r, c, self.clues[(r, c)])]
            self.dsu.init_root_stats(
                (r, c, side),
                rat=subcell_rational(st, side),
                pi_q=subcell_pi_count(st, side),
                is_near=(side == NEAR),
                is_far=(side == FAR),
                cell=(r, c),
                clue=clue,
                init_open=open_init,
            )

        # Process edges to ASSIGNED neighbors:
        # - decrement that neighbor's sub-cell root's open_edges (was +1, now 0).
        # - union the two sub-cells.
        for side in subcells(st):
            edges_touched = subcell_edges(st, side)
            for edge in edges_touched:
                dr, dc = EDGE_DELTA[edge]
                nr, nc = r + dr, c + dc
                if not (0 <= nr < n and 0 <= nc < n):
                    continue
                if not self._is_assigned(nr, nc, idx):
                    continue
                nst = states[nr][nc]
                nedge = OPPOSITE_EDGE[edge]
                for nside in subcells(nst):
                    if nedge in subcell_edges(nst, nside):
                        # The neighbor sub-cell's edge into (r,c) was open;
                        # now it's closed (cell (r,c) just got assigned).
                        # The (r,c) sub-cell's edge into the neighbor was
                        # already counted as NOT open in init_open above
                        # (since neighbor was assigned).
                        self.dsu.add_open_edges((nr, nc, nside), -1)
                        # Union.
                        self.dsu.union((r, c, side), (nr, nc, nside))

        # No-dangling: if this cell is arced and its NEAR / FAR are in the same DSU
        # component, immediate fail.
        if st != NO_ARC:
            if self.dsu.find((r, c, NEAR)) == self.dsu.find((r, c, FAR)):
                return False
        return True

    def _after_add_check(self, states, r, c, idx):
        """Run the cheap interior-vertex check, then validate any newly-closed
        regions (open_edges hit zero) and apply clue-area-bound pruning to
        open regions containing clue cells.
        """
        n = self.n

        # Interior vertex no-dangle check (vertex (r, c) finalized).
        if r >= 1 and c >= 1 and r <= n - 1 and c <= n - 1:
            cnt = _endpoint_count_at_vertex(states, r, c, n)
            if cnt == 1:
                return False

        # Examine the regions touching the newly-added cell. They are the
        # roots of the newly-added subcells (after unions).
        st = states[r][c]
        roots_touched = set()
        for side in subcells(st):
            roots_touched.add(self.dsu.find((r, c, side)))

        for root in roots_touched:
            # Open-region prunes:
            # No-dangle (between near and far cell sets, via merged regions).
            if self.dsu.near_cells[root] & self.dsu.far_cells[root]:
                return False
            # Clue area lower-bound.
            for (cr, cc, label) in self.dsu.clues[root]:
                if self.dsu.rat[root] > label:
                    return False
            # Closed region?
            if self.dsu.open_edges[root] == 0:
                # Integer area.
                if self.dsu.pi_q[root] != 0:
                    return False
                area = self.dsu.rat[root]
                # If clue, check pieces*area = label.
                clues = self.dsu.clues[root]
                if clues:
                    members = self.dsu.members[root]
                    segs = collect_region_segments(states, members, n)
                    pieces = count_smooth_pieces(segs)
                    for (cr, cc, label) in clues:
                        if area * pieces != label:
                            return False

        # Some neighbor regions might have been closed by the open-edge
        # decrement (their open_edges hit 0) even if they don't contain
        # subcells of the newly-added cell. Check them too.
        # Find neighbor sub-cells we just decremented:
        for edge, (dr, dc) in EDGE_DELTA.items():
            nr, nc = r + dr, c + dc
            if not (0 <= nr < n and 0 <= nc < n):
                continue
            if not self._is_assigned(nr, nc, idx):
                continue
            nst = states[nr][nc]
            for nside in subcells(nst):
                root = self.dsu.find((nr, nc, nside))
                if root in roots_touched:
                    continue  # already validated above
                roots_touched.add(root)
                # Same checks.
                if self.dsu.near_cells[root] & self.dsu.far_cells[root]:
                    return False
                for (cr, cc, label) in self.dsu.clues[root]:
                    if self.dsu.rat[root] > label:
                        return False
                if self.dsu.open_edges[root] == 0:
                    if self.dsu.pi_q[root] != 0:
                        return False
                    area = self.dsu.rat[root]
                    clues = self.dsu.clues[root]
                    if clues:
                        members = self.dsu.members[root]
                        segs = collect_region_segments(states, members, n)
                        pieces = count_smooth_pieces(segs)
                        for (cr, cc, label) in clues:
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
