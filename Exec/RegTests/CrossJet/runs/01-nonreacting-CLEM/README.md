# 01 — non-reacting CLEM
Input: `inputs-case1-nonreacting-CLEM.inp`
Config: `pelec.do_les=0`, `pelec.Cs=0.16`, `clem.do_react=0`, `clem.do_stir=1`
Data: 19 plotfiles, t = 0 -> 0.18 ms. No checkpoints retained.
Note: with `do_les=0`, CLEM's stirring uses Smagorinsky nu_t (nonzero at the
wall) while NO SGS momentum term is added. This is the cleanest momentum
treatment of the four runs.
