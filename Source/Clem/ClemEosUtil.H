#ifndef CLEMEOSUTIL_H
#define CLEMEOSUTIL_H

#include "PelePhysics.H"

#include "ClemIndex.H"

namespace clem::eos_util {

// Fernando-Clem: bias-free (rho, e, Y) -> T inversion, ported from the old
// implementation. Plain REY2T extrapolates linearly below its 300 K bracket
// (unburned gas sits below it) and that bias, converted back to energy every
// step, acts as a spurious heating ratchet. A short Newton polish with
// forward RTY2E/RTY2Cv evaluations removes it.
AMREX_FORCE_INLINE
void
REY2T_consistent(
  const amrex::Real rho,
  const amrex::Real e_spec,
  const amrex::Real* Y,
  amrex::Real& T)
{
  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real Yloc[NUM_SPECIES];
  for (int sp = 0; sp < NUM_SPECIES; sp++) {
    Yloc[sp] = Y[sp];
  }

  eos.REY2T(rho, e_spec, Yloc, T); // bracketed / extrapolated first estimate

  for (int it = 0; it < 3; ++it) {
    amrex::Real eT = 0.0;
    amrex::Real cv = 0.0;
    eos.RTY2E(rho, T, Yloc, eT);
    eos.RTY2Cv(rho, T, Yloc, cv);
    const amrex::Real dT = (e_spec - eT) / cv;
    T += dT;
    if (std::abs(dT) < 1e-12 * T) {
      break;
    }
  }
}

} // namespace clem::eos_util
#endif
