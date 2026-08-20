# 02 — reacting CLEM + constant Smagorinsky  [BROKEN — keep for the record]
Input: `inputs-case1-reacting-CLEM.inp`
Config: `do_les=1`, `les_model=0` (default), `Cs=0.16`, `PrT=1.0`, `clem.do_react=1`
Data: 67 plotfiles t = 0 -> 0.66 ms, 14 checkpoints, `datlog`.
Result: undamped Smagorinsky gives nu_t ~ 50x molecular in the first wall cell.
Boundary layer 7x too thick, displacement thickness 21x too large, upstream
mean pressure +19.3%, oblique wave lattice reflecting between the walls, wall
cell at M = 0.10, flame walks to the inlet at 100-140 m/s.
Do not use for physics conclusions.
