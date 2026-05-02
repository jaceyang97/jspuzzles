"""Entry point: solve a puzzle and report the answer."""

import argparse
import time
from puzzle import (
    N, GREEN, CLUES,
    EXAMPLE_N, EXAMPLE_GREEN, EXAMPLE_CLUES, EXAMPLE_ANSWER,
)
from solver import Solver
from score import evaluate
from geometry import NO_ARC


def fmt_states(states):
    sym = {0: '.', 1: '\\', 2: '/', 3: '/', 4: '\\'}
    sym_full = {0: 'none', 1: 'NW', 2: 'NE', 3: 'SW', 4: 'SE'}
    out = []
    for row in states:
        out.append(' '.join(f'{sym_full[s]:>4}' for s in row))
    return '\n'.join(out)


def run(n, green, clues, *, expected=None, max_steps=None, verbose=True, max_solutions=1, order='row_major'):
    print(f'Puzzle {n}x{n}: {len(clues)} clues, {len(green)} green cells. order={order}')
    s = Solver(n, green, clues, verbose=verbose)
    t0 = time.time()
    try:
        s.solve(max_solutions=max_solutions, max_steps=max_steps, order_kind=order)
    except KeyboardInterrupt:
        print('Interrupted')
    elapsed = time.time() - t0
    print(f'Steps: {s.steps:,}  Elapsed: {elapsed:.2f}s  Solutions: {len(s.solutions)}  max_idx={s.max_idx}/{n*n}')
    for k, (states, info) in enumerate(s.solutions):
        print(f'\n--- Solution {k} ---')
        print(fmt_states(states))
        res = evaluate(states, clues, strict=False)
        print('\nCell scores:')
        for row in res['cell_score']:
            print(' ', ' '.join(f'{x:>4}' for x in row))
        print(f"row_sums={res['row_sums']}")
        print(f"col_sums={res['col_sums']}")
        print(f"ANSWER = {res['answer']}")
        if expected is not None:
            print(f"expected = {expected}  ", '✓' if res['answer'] == expected else '✗')
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--example', action='store_true', help='Run the 4x4 example')
    ap.add_argument('--max-steps', type=int, default=None)
    ap.add_argument('--quiet', action='store_true')
    ap.add_argument('--order', default='row_major', choices=['row_major', 'boundary_first', 'clues_first'])
    args = ap.parse_args()

    if args.example:
        run(EXAMPLE_N, EXAMPLE_GREEN, EXAMPLE_CLUES,
            expected=EXAMPLE_ANSWER, max_steps=args.max_steps,
            verbose=not args.quiet, order=args.order)
    else:
        run(N, GREEN, CLUES, max_steps=args.max_steps, verbose=not args.quiet, order=args.order)


if __name__ == '__main__':
    main()
