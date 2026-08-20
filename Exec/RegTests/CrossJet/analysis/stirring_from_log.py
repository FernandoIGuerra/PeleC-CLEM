"""Pull the per-timestep CLEM stirring counters out of a PeleC run log.

The stirring diagnostic is printed once per step by ClemAlgorithm.cpp:

    CLEM stirring: maps applied = 112 over 307200 lines (max = 1/line),
    eddies rejected as too small = 76648 (99.85% of those drawn),
    lines truncated by the substep cap = 0

It carries no step number of its own, so each one is tied to the most recent

    [Level 0 step N] ADVANCE at time T with dt = DT

which opens the same step block. dt is what makes the closures comparable:
the two runs take different step sizes, so "maps per step" and "maps per
second of physical time" rank them differently - see the rate column.

    python3 stirring_from_log.py run.log out.csv

Safe to run against a log that is still being written; the trailing partial
block is just dropped.
"""
import csv
import re
import sys

ADVANCE = re.compile(
    r"^\[Level 0 step (\d+)\] ADVANCE at time (\S+) with dt = (\S+)")
STIR = re.compile(
    r"^CLEM stirring: maps applied = (\d+) over (\d+) lines "
    r"\(max = (\d+)/line\), eddies rejected as too small = (\d+) "
    r"\(([-+0-9.eE]+)% of those drawn\), "
    r"lines truncated by the substep cap = (\d+)")

COLS = ["step", "time", "dt", "maps", "lines", "max_per_line",
        "rejected", "rejected_pct", "truncated", "drawn", "maps_per_sec"]


def parse(path):
    rows, step, t, dt = [], None, None, None
    with open(path, errors="replace") as f:
        for line in f:
            if line.startswith("[Level 0 step "):
                m = ADVANCE.match(line)
                if m:
                    step, t, dt = int(m[1]), float(m[2]), float(m[3])
            elif line.startswith("CLEM stirring: maps applied"):
                m = STIR.match(line)
                if m and step is not None:
                    maps, lines_, mx = int(m[1]), int(m[2]), int(m[3])
                    rej, pct, trunc = int(m[4]), float(m[5]), int(m[6])
                    rows.append([step, t, dt, maps, lines_, mx, rej, pct,
                                 trunc, maps + rej,
                                 maps / dt if dt else 0.0])
                    step = None      # one stirring line per ADVANCE
    return rows


def main():
    src, dst = sys.argv[1], sys.argv[2]
    rows = parse(src)
    if not rows:
        sys.exit(f"no stirring records found in {src}")
    with open(dst, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(COLS)
        w.writerows(rows)

    maps = [r[3] for r in rows]
    rate = [r[10] for r in rows]
    nz = [m for m in maps if m > 0]
    tot_rej, tot_drawn = sum(r[6] for r in rows), sum(r[9] for r in rows)
    print(f"{src}")
    print(f"  steps with a stirring record : {len(rows)}")
    print(f"  physical time covered        : {rows[0][1]:.6e} -> "
          f"{rows[-1][1]:.6e} s")
    print(f"  maps applied / step          : mean {sum(maps)/len(maps):.2f}, "
          f"max {max(maps)}, total {sum(maps)}")
    print(f"  steps with zero maps         : {len(maps) - len(nz)} "
          f"({100.0 * (len(maps) - len(nz)) / len(maps):.1f}%)")
    if nz:
        print(f"  mean over non-zero steps     : {sum(nz)/len(nz):.2f}")
    print(f"  maps / second (physical)     : mean "
          f"{sum(rate)/len(rate):.4e}, max {max(rate):.4e}")
    print(f"  eddies rejected as too small : {tot_rej} of {tot_drawn} drawn "
          f"({100.0 * tot_rej / tot_drawn if tot_drawn else 0:.3f}%)")
    print(f"  wrote {dst}")


if __name__ == "__main__":
    main()
