# CrossJet — H2 jet in supersonic crossflow, Case 1

Jin, Cai, Hong, Zhang & Liang, *Combustion and Flame* **268** (2024) 113627.
Domain 12.5 x 1.5 x 0.5 cm, grid 400 x 48 x 16 (dx = dy = dz = 0.3125 mm),
single level, BurkeDryer / Fuego / Simple transport.

Last updated after the run-set below was completed and the directory
reorganised. **No jobs are running.**

## Directory layout

    runs/01-nonreacting-CLEM/            plt00000-14992
    runs/02-reacting-CLEM-smagorinsky/   plt_React*, chk*, datlog
    runs/03-reacting-CLEM-WALE-wallmodel/    plt_WM*, chk_wm*
    runs/04-reacting-CLEM-WALE-nowallmodel/  plt_WALE*, chk_wale*
    analysis/    figures, plotting scripts, digitised paper data
    paper/       the paper PDF and case data
    _trash/      plt_clem_init (t = 0 CLEM startup dump, not data)

Each run directory carries a copy of the input it was produced with.
The **no-model reference** lives in a different tree:
`/home/pc08/PeleC/PeleC/Exec/RegTests/CrossJet/plt_React*` (t = 0 -> 1.0 ms).

## The runs

| # | configuration | t reached | ignition | note |
|---|---|---|---|---|
| ref | pristine PeleC, no model, `do_react=1` | 1.00 ms | **0.63-0.70 ms, x = 12.48 cm (outlet, wall cell)** | matches the paper's mechanism |
| 01 | CLEM, `do_les=0`, `clem.do_react=0` | 0.18 ms | n/a (non-reacting) | first run |
| 02 | CLEM + **Smagorinsky** (`les_model=0`, Cs=0.16) | 0.66 ms | yes | **BROKEN** - see below |
| 03 | CLEM + **WALE** + **wall model** | 0.62 ms | 0.36-0.42 ms at the outlet, then walks upstream past the jet | inlet waves fixed, flame does not anchor |
| 04 | CLEM + **WALE**, no wall model | 0.94 ms | **NONE** - maxT flat ~1215 K, Y(H2O) ~2e-5 | under-mixed, see below |

## What was established

1. **Run 02's inlet wave lattice was constant-Smagorinsky, not CLEM.**
   Undamped `les_model=0` puts nu_t ~ 50x molecular in the first cell at a
   no-slip wall (nu_t = Cs^2 Delta^2 |S| ~ 74 vs nu ~ 1.2 cm^2/s). Boundary
   layer 7x too thick, displacement thickness 21x too large, upstream mean
   pressure +19.3% vs +1.0%. That launched oblique waves at the Mach angle
   which reflected between the walls. CLEM cannot cause this: it writes only
   rho Y_k and rho e, never momentum.
   Note upstream PeleC *aborts* for NUM_SPECIES > 2 with LES; that guard is
   commented out in this fork, which is why the configuration was reachable.

2. **The upstream flame follows from a subsonic near-wall channel.**
   A flame cannot cross a M = 1.5 core; it crawls up a subsonic wall layer.
   Run 02 had M = 0.10 in the wall cell at x = 3 cm at t = 0.18 ms (reference:
   M = 1.32, no subsonic layer anywhere). The channel forms BEFORE ignition,
   so it is cause, not consequence.

3. **WALE fixes the pressure lattice** (run 03) - the wave field is gone and
   the upstream pressure returns to near-ambient.

4. **A WALE bug was found and fixed** in `Source/LES.H`
   (`get_wale_sfs_stresses`): the `mut` denominator had no zero guard, so a
   uniform region gives 0/0 = NaN, which propagates into the momentum and
   energy SFS fluxes. Confirmed with `amrex.fpe_trap_invalid=1` (traps at
   LES.H:496, no CLEM in the stack). **Still present in upstream PeleC**
   (ac17bd0, 2026-07-27, Source/LES.H:494); their CI misses it because the
   NUM_SPECIES > 2 abort confines WALE to TG/HIT, where every cell has a
   gradient.

5. **Run 04 does not ignite at all**, to 0.94 ms - 0.3 ms past the reference's
   ignition, with no precursor (reference at 0.58 ms already had
   Y(H2O) = 0.03). The likely cause is a real coupling gap: WALE returns
   nu_t = 0 identically in pure shear, CLEM takes nu_t from the active LES
   model, so `kolmogorovScale -> 0`, `eddyTime -> never`, and **no triplet maps
   fire in the near-wall cells** - exactly where the reference ignites
   (y = 0.02 cm). Run 03's wall model masked this by thickening the layer.
   Over-mixing (run 02) has been traded for under-mixing (run 04).

## Open items

* **Next run**: `pelec.do_les = 0`, no wall model, Maxwell CLEM unchanged.
  Because `les.model = do_les ? les_model : 0` (PeleC.cpp:1041), this pins
  CLEM's stirring to Smagorinsky with Cs = 0.16 - **nonzero at the wall**, so
  near-wall stirring returns - while adding no SGS momentum term at all. It is
  also the configuration closest to the paper, which uses **no turbulence
  model** ("without incorporating a turbulent model").
* **Wall model caveat** (`Source/WallModel.H`): the Spalding law amplifies
  tau_w 4.6x at these conditions because y1+ ~ 32 never sees the viscous
  sublayer, while the inflow is a slug with no boundary layer. Solver itself is
  unit-tested (root residual ~1e-16, scale = 1.000 in the laminar limit, log-law
  asymptote to 0.3%).
* **nu_t floor**: van Driest `nu_t = max(nu_t_WALE, kappa*y_w*u_tau*
  [1-exp(-y+/26)]^2)` from the wall model's u_tau would close the run-04 gap.
  Not implemented.
* **Resolution**: this grid is the paper's *coarse* L2 (y+ = 25 there, ~32
  measured here) with no AMR; the paper's production is L3 (finest 0.0781 mm,
  y+ = 12) and it rejects L2 as under-resolved. CLEM only runs on level 0
  (`clem_active = (level == 0 && ...)`, Advance.cpp:77), so AMR is not
  available with CLEM on.
* **NLEM**: raising `CLEM_NLEM` requires `-DNUM_LEM` to match in PelePhysics
  (`ReactorClemConfig.H`); `Make.PeleC` propagates only `CLEM_NLEM`, so a
  mismatch is a compile error.

## Reproducing the analysis

`analysis/` holds the plotting scripts. They contain a standalone single-level
AMReX plotfile reader, so no yt is needed. Use the only python with matplotlib:

    /home/pc08/miniconda3/envs/felics/bin/python3

Paths inside those scripts predate this reorganisation and point at the old
flat layout - update them to `runs/<case>/` before rerunning.
