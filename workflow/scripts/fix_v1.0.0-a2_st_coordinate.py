import sys
import polars as pl

from typing import Any
from itertools import islice
from collections import deque


def sliding_window(iterable, n):
    "Collect data into overlapping fixed-length chunks or blocks."
    iterator = iter(iterable)
    window = deque(islice(iterator, n - 1), maxlen=n)
    for x in iterator:
        window.append(x)
        yield tuple(window)


def main():
    infile = sys.argv[1]
    # Edge case where one interval is null interval. Need to double check.
    df = pl.read_csv(infile, separator="\t", has_header=True).filter(
        pl.col("chromStart") != pl.col("chromEnd")
    )
    fixed_itvs = set()
    print(*df.columns, sep="\t", file=sys.stdout)

    for r1, r2 in sliding_window(df.iter_rows(), 2):
        if r1 in fixed_itvs:
            continue

        # Last interval start will be +1 in v1.0.0-a2
        if r1[0] == r2[0] and r1[2] != r2[1]:
            # Add original interval so won't redo.
            fixed_itvs.add(r2)
            r2 = (r2[0], r2[1] - 1, *r2[2:])
            print(*r1, sep="\t", file=sys.stdout)
            print(*r2, sep="\t", file=sys.stdout)
        else:
            print(*r1, sep="\t", file=sys.stdout)

    final_row = df.row(-1)
    if final_row not in fixed_itvs:
        print(*final_row, sep="\t", file=sys.stdout)


if __name__ == "__main__":
    raise SystemExit(main())
