#ifndef CLEM_NUMERICS_CONFIG_H
#define CLEM_NUMERICS_CONFIG_H

// ============================================================
// LEM thermal diffusion scheme — set at compile time:
//   make CLEM_DIFF_SCHEME=<int>        (GNUmake)
//   cmake -DCLEM_DIFF_SCHEME=<int>     (CMake, if wired)
//
//   0  CLEM_DIFF_SCHEME_O2_EULER   : 2nd-order space, explicit Euler (baseline)
//   1  CLEM_DIFF_SCHEME_O4_EULER   : 4th-order space, explicit Euler
//   2  CLEM_DIFF_SCHEME_O2_RK2     : 2nd-order space, Heun's RK2
//   3  CLEM_DIFF_SCHEME_O2_SUBSTEP : 2nd-order space, CFL-based sub-stepping (default)
//                                     Computes dt_cfl = dm^2/(2*beta_max), takes
//                                     N_sub = ceil(dt/dt_cfl) explicit Euler steps,
//                                     and passes the averaged rate to the reactor.
//                                     Unconditionally stable for any dt.
//
// Stability of single-step schemes (r = alpha_eff*dt/dm^2 must be <= 0.5):
//   - O2_EULER   : r <= 0.5
//   - O4_EULER   : r <= 3/11 ≈ 0.27  (wider stencil raises max eigenvalue)
//   - O2_RK2     : r <= 0.5           (same stability as Euler for diffusion)
//   - O2_SUBSTEP : unconditionally stable via automatic sub-stepping
// ============================================================
#define CLEM_DIFF_SCHEME_O2_EULER   0
#define CLEM_DIFF_SCHEME_O4_EULER   1
#define CLEM_DIFF_SCHEME_O2_RK2     2
#define CLEM_DIFF_SCHEME_O2_SUBSTEP 3

#ifndef CLEM_DIFF_SCHEME
#define CLEM_DIFF_SCHEME CLEM_DIFF_SCHEME_O2_SUBSTEP
#endif

#endif /* CLEM_NUMERICS_CONFIG_H */
