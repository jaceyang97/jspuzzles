# Arc-drawing puzzle solver

Solver framework for the "draw 90° arcs in some white cells" puzzle.

## Modules

| File          | Purpose                                                                   |
|---------------|---------------------------------------------------------------------------|
| `puzzle.py`   | Puzzle data: grid size, green cells, numbered clues. Edit to suit.        |
| `geometry.py` | Arc encoding (5 states per cell), sub-cell decomposition, edge adjacency. |
| `regions.py`  | Connected-components finder over sub-cells.                               |
| `boundary.py` | Boundary-segment collector and smooth-piece counter for a region.         |
| `score.py`    | Region score (= smooth-pieces × area), full puzzle evaluation.            |
| `solver.py`   | Backtracking solver (incremental DSU + rollback).                         |
| `main.py`     | CLI entry point.                                                          |

## Status

- 4x4 example: **solved correctly**, answer = **18,928** in 3,001 steps (<1s).
- 9x9 puzzle: search reaches ~70% depth (max_idx ≈ 60/81) within seconds but
  the remaining branches are too many for plain backtracking. With 5^63 ≈ 7×10^43
  raw configurations and the current ~10^5x average pruning factor, we still
  have ~10^38 nodes to explore — infeasible.

## Constraints checked

1. **Greens** never have arcs.
2. **No-dangling**: an arced cell's two parts must end up in different regions.
3. **Integer area**: each region's π/4 contributions must cancel
   (= equal NEAR and FAR sub-cell counts in the region).
4. **Clue match**: for each numbered cell, the region containing its larger
   side has score = label.
5. **Interior-vertex no-dangle**: at any interior lattice vertex, the count
   of arc endpoints meeting there must not equal 1 (Y-junctions of degree 3
   are fine — they correspond to three regions meeting at a point).

## What's missing for the 9x9

Brute-force backtracking is too slow. Possible directions:

- **Stronger constraint propagation**: domain reduction across lattice
  vertices, region-area divisibility, smooth-piece upper bounds from
  clue label / area divisor pairs.
- **Constraint solver**: encode the puzzle into MiniZinc / CPLEX / Z3 with
  region-connectivity auxiliary variables.
- **Region-first search**: enumerate valid region shapes around each clue
  cell (factorizations of the clue label into pieces × area), then verify
  global consistency.

## Usage

```sh
python3 main.py --example       # 4x4 example, expect 18,928
python3 main.py                 # 9x9 puzzle (will run a long time)
python3 main.py --max-steps 10000000 --order boundary_first
```

The `--order` flag selects variable ordering: `row_major` (default),
`boundary_first` (outer ring first), or `clues_first` (numbered cells first).
