# 04 — reacting CLEM + WALE, NO wall model  (isolation run)
Input: `inputs-case1-reacting-CLEM-WALE.inp`
Config: identical to run 03 except `pelec.do_wall_model = 0`
Data: 95 plotfiles t = 0 -> 0.94 ms, 8 checkpoints.
Result: **no ignition at all** — max T flat at ~1215 K, max Y(H2O) ~2e-5, at
0.94 ms. The reference had ignited by 0.70 ms and reached 2336 K by 1.0 ms.
Interpretation: WALE gives nu_t = 0 identically in pure shear, CLEM takes its
nu_t from the active LES model, so no triplet maps fire in the near-wall cells
— which is exactly where the reference ignites (y = 0.02 cm). Under-mixed.

CAVEAT: two 64-rank jobs were briefly running on this input at the same time
and wrote the same filenames. The surviving series has 95 files at exactly
1e-5 s spacing with no duplicate times, so it looks like one clean trajectory,
but individual snapshots cannot be *proven* to belong to a single realization
(CLEM stirring draws amrex::Random()). The no-ignition conclusion is
unaffected — neither realization ignited.
