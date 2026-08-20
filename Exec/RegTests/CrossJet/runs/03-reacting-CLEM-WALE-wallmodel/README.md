# 03 — reacting CLEM + WALE + wall model
Input: `inputs-case1-reacting-CLEM-WM.inp`
Config: `les_model=2` (WALE, `Cw=0.325`), `PrT=1e6`, `do_wall_model=1`
Data: 63 plotfiles t = 0 -> 0.62 ms, 13 checkpoints.
Resume with `amr.restart=chk_wm52945` (t = 6.0e-4 s).
Result: inlet pressure lattice GONE. Ignites 0.36-0.42 ms at the outlet
(x ~ 12.3 cm) — same location as the reference and the paper — but too early
(reference: 0.63-0.70 ms) and the flame does not anchor at the injector; it
propagates past it to x ~ 2.5 cm. Wall cell still subsonic (M = 0.06),
crossing M = 1 at x ~ 1.7 cm.
