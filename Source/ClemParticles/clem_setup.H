#ifndef CLEM_SETUP_H
#define CLEM_SETUP_H

#include "ClemParticles/clem_numerics_config.H"

// ============================================================
// CLEM user setup — edit this file to configure the simulation.
// All choices are compile-time constants (constexpr).
// ============================================================

namespace clem {
namespace setup {

    /**
     * Thermal diffusion scheme for the LEM periodic ring.
     * See clem_numerics_config.H for the full description of each option.
     *
     *   DiffusionScheme::O2_Euler   — 2nd-order space / explicit Euler
     *   DiffusionScheme::O4_Euler   — 4th-order space / explicit Euler
     *   DiffusionScheme::O2_RK2     — 2nd-order space / Heun's RK2
     *   DiffusionScheme::O2_Substep — 2nd-order space / CFL sub-stepping  ← recommended
     */
    constexpr DiffusionScheme diffusion_scheme = DiffusionScheme::O2_Substep;

} // namespace setup
} // namespace clem

#endif /* CLEM_SETUP_H */
