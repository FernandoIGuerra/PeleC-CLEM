# CLEM 1D premixed flame — accuracy test

A 1D freely-propagating premixed flame used to evaluate how accurately the
embedded LES-LEM (CLEM particle) model predicts:

1. **displacement / laminar flame speed** `S_L`
2. **flame thickness** `delta_th`
3. **production (burning) rate** of the flame

## Why a 1D flame is the right test

Each LES cell carries an LEM line (`NUM_LEM` sub-points) on which diffusion +
reaction are solved; transport across the domain happens through advection and
splicing of the particles. With **no turbulent stirring (no triplet maps)** an
LEM line reduces *exactly* to the 1D laminar reaction–diffusion flame equation.
So in the laminar limit the CLEM model must reproduce a detailed-chemistry
laminar flame — same speed, thickness and reaction-rate profile. This isolates
the diffusion discretization (`clem_diffusion.cpp`), the reaction coupling, and
the particle advection/splicing chain, with no turbulence confound.

## Configuration

Stabilized (stationary) flame, quasi-1D:

- transverse periodic, **inflow (`Hard`)** of cold reactants at `x-lo`,
  **outflow (`FOExtrap`)** at `x-hi`
- mixture: **phi = 0.4 H2/air, LiDryer, p = 1 atm, Tu = 298 K**
- reference: `LiDryer_H2_p1_phi0_4000tu0300.dat` (PREMIX/Cantera), `S_L ≈ 22.8 cm/s`

The **whole domain is initialized from the reference flame profile** (PeleC PMF
reader): the flame starts at its correct thickness, structure and velocity field
(cold-end velocity = S_L is also the Hard inflow), so only fast acoustic
transients (≈ L/c, microseconds) need to settle — *not* a full convective
flow-through. `prob.standoff` is auto-set so the reference front lands at
`prob.x_flame`; override it to move the flame.

This makes the test a **maintenance / stationarity test**: does CLEM *hold* the
correct S_L and δ_th, or does numerical diffusion thicken the flame over a few
flame times (δ²/α ≈ 2.7e-3 s)?  Watch δ_th over successive plotfiles — drift
upward = numerical broadening.

> Note: brute-force convergence from a generic (e.g. tanh) IC is impractical
> here because the flame Mach number is tiny (Ma ≈ S_L/c ≈ 6.5e-4): PeleC takes
> an acoustic timestep, so steady state needs ~millions of steps. Starting from
> the true profile sidesteps that.

## Measurables (computed by `postprocess.py`)

- `S_L = mean(rho*u)/rho_u` — mass flux is constant across a steady flame, so
  this is frame-independent and robust.
- `delta_th = (T_b - T_u)/max|dT/dx|`
- burning rate `= rho_u * S_L * (Y_H2,u - Y_H2,b)` = ∫ H2 consumption rate dx.

```
make -j                       # or cmake build
./PeleC2d.gnu.MPI.ex inputs.inp
python postprocess.py CurrentResults/pltFlame_XXXXX
```

## Turning this into an accuracy study

1. **Reference truth** is the in-repo `.dat` (matched mixture/p/Tu). Optionally
   regenerate with Cantera for other conditions.
2. **Apples-to-apples second reference:** run the same mixture/mesh with native
   PeleC (the `PMF` RegTest, fully resolved). CLEM−PMF difference = model error
   separated from mechanism/EOS error.
3. **Convergence sweep (headline result):** refine intra-line resolution
   `NUM_LEM` (currently 9, in `Source/ClemParticles/clem_index_definition.H`)
   and/or LES `dx` (`amr.n_cell`). Show `S_L`, `delta_th`, burning rate
   converging to the reference; report error vs. `NUM_LEM` and `dx`.
4. **Scheme comparison:** repeat the sweep for `O2_Euler` / `O4_Euler` /
   `O2_RK2` (set in `Source/ClemParticles/clem_setup.H`) to show observed order
   of accuracy and which converges fastest — direct validation of the diffusion
   discretization work.

## Cross-check: freely propagating (transient)

Set `u_in = 0` (or a small value), ignite near `x-lo`, and track the front
position over time; `S_L = dx_front/dt − u_unburned`. Independent check of the
stationary-flame `S_L`.
