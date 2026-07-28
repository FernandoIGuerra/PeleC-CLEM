#ifndef CLEMREACTION_H
#define CLEMREACTION_H

#include <array>

#include "PelePhysics.H"

#include "ClemManager.H"
#include "ClemScratch.H"

namespace pele::physics::reactions {
class ReactorCvode;
} // namespace pele::physics::reactions

namespace clem {

// Fernando-Clem: subgrid chemistry along the LEM line (Maxwell [M 3.45-3.46]:
// DT/Dt = -(gamma-1)/gamma * Q*omega/rho, DY_k/Dt = omega_k/rho), integrated
// per element by the PERSISTENT, per-rank CVODE reactor
// (pele::physics::reactions::ReactorCvode::reactClemLine - see
// PelePhysics/Source/Reactions/ReactorCvode.H/.cpp). All config::n_lem
// elements of one line are reacted as ONE batched CVODE system
// (icell = lem_index), so the whole line costs one CVodeReInit+CVode call,
// not n_lem separate ones, and the CVODE context itself (matrix, linear/
// nonlinear solver, user data) is built ONCE per rank on first use (reactor(),
// a Meyer's singleton shared by BOTH entry points below) and reused for the
// rest of the run.
//
// SPLITTING - SEQUENTIAL (Lie), not coupled: reaction always runs AFTER
// diffusion has updated the line, with the reactor's external source (rhs)
// forced to ZERO, never fed from diffusion's rate. This is chosen over
// coupling through CVODE's rhs because the particle state is PRIMITIVE (T,
// Y_k, not rho*Y_k/rho*e) - building a properly CONSERVATIVE external source
// from a primitive line would itself be an approximation.
//
// GRANULARITY - two entry points, chosen by whether diffusion is bounding the
// step (Maxwell [M §3.3.3] steps 14-20: Delta t_diff comes from the diffusion
// CFL; diffusion and reaction are integrated TOGETHER over that substep, not
// as two independent full-Delta_t stages):
//   - reactSubstep(): do_diffusion=1. Called from DiffusionOperator::diffuse's
//     OWN CFL-limited substep loop, immediately after each advanceLine, on
//     the SAME dt_sub. Operates directly on the LineScratch diffusion already
//     gathered - no redundant gather/scatter round-trip to the particles
//     between the two operators.
//   - react(): do_diffusion=0. There is no diffusion CFL to bound a substep -
//     Maxwell's Delta t_diff degenerates to the full LES Delta t - so this
//     gathers/reacts/scatters every line once, over the whole dt, directly
//     against the particle AoS.
//
// ISOBARIC, matching the diffusion stage's closure and Maxwell's constant-
// pressure LEM array ([M §3.2.2]): the reactor runs in HP mode
// (ReactorTypes::h_reactor_type) - conservatives[][NUM_SPECIES] carries
// rho*h (specific enthalpy density), not rho*e, and CVODE's RHY2T inversion
// recovers T. Density is NOT frozen: after reaction, (rho, e) are
// RECONSTRUCTED from the new (p_held, T, Y) via PYT2RE - the exact same
// isobaric closure DiffusionOperator::advanceLine uses - so pressure is left
// unchanged (still the reference p_1 the next isentropic LES->SGS stage
// compresses from) and the line's volume breathes with the heat release.
//
// NUM_LEM (PelePhysics, ReactorClemConfig.H) sizes the persistent reactor's
// CVODE vector at build time and MUST equal clem::config::n_lem - enforced
// by a static_assert in ClemReaction.cpp, where both are visible together.

// Fernando-Clem: per-rank reusable scratch for the reaction operator, exact
// same pattern as ClemScratch.H's LineScratch (grown/owned once per rank,
// reused unchanged for every line of every step afterwards) - here the size
// is compile-time fixed (config::n_lem, no dynamic resize needed), so plain
// std::array on a thread_local static instance is enough, no heap at all.
struct ReactionScratch
{
  static ReactionScratch& get()
  {
    static thread_local ReactionScratch arena;
    return arena;
  }

  std::array<int, config::n_lem> pidx;
  std::array<std::array<amrex::Real, NUM_SPECIES + 1>, config::n_lem>
    primitives;
  std::array<std::array<amrex::Real, NUM_SPECIES + 1>, config::n_lem>
    conservatives;
  std::array<std::array<amrex::Real, NUM_SPECIES + 1>, config::n_lem> rhs;
  std::array<amrex::Real, config::n_lem> density;
  std::array<amrex::Real, config::n_lem> press_held;
  std::array<amrex::Real, config::n_lem> mass;
};

class ReactionOperator
{
public:
  // Fernando-Clem: diagnostics of one reaction pass (MPI-reduced)
  struct Diagnostics
  {
    amrex::Long n_lines = 0;     // lines reacted
    amrex::Long n_skipped = 0;   // lines skipped (empty cell / placeholder)
    amrex::Real max_neg_Y = 0.0; // worst Y clipped by clean_init_massfrac
    // Fernando-Clem: worst per-line substep-halving retry count from
    // ReactorCvode::reactClemLine's while-loop (0 = every line reached the
    // full dt in one CVode() call). Nonzero is not an error - the loop
    // guarantees full integration or a loud abort - but a persistently high
    // count flags stiff lines worth a closer look (mechanism/dt/tolerances)
    int max_retries = 0;
    amrex::Long n_retried_lines = 0; // lines that needed >=1 retry
  };

  // Fernando-Clem: react every line on the level by dt, directly against the
  // particle AoS. Standalone (do_diffusion=0) path only - see the file
  // comment's GRANULARITY note; when diffusion is on, reaction is co-stepped
  // via reactSubstep() instead, called from inside diffuse()'s substep loop.
  static Diagnostics react(ClemManager& mgr, int lev, amrex::Real dt);

  // Fernando-Clem: react ONE already-gathered line over dt_sub - the
  // diffusion-CFL substep. Called from DiffusionOperator::diffuse's own
  // substep loop, immediately after each advanceLine, on the SAME dt_sub
  // (Maxwell [M §3.3.3] steps 14-20: diffuse+react together over the
  // CFL-limited substep, repeated until the full LES dt is covered).
  // Operates directly on lscr's (T, rho, eint, Y) - the same arena diffusion
  // gathered/will scatter, no round-trip through the particles in between.
  // Updates diag's max-reduced fields (max_neg_Y, max_retries) directly;
  // does NOT touch n_lines/n_retried_lines (those are per-CELL counts - the
  // caller increments them once after the whole substep loop for a line
  // finishes). Returns the retry count from this substep's reactClemLine
  // call (0 = clean).
  static int
  reactSubstep(LineScratch& lscr, amrex::Real dt_sub, Diagnostics& diag);

private:
  // Fernando-Clem: the persistent, per-rank CVODE reactor - a Meyer's
  // singleton so react() and reactSubstep() share the SAME context (built
  // once, lazily, on whichever entry point is called first) rather than each
  // constructing its own.
  static pele::physics::reactions::ReactorCvode& reactor();

  // Fernando-Clem: gather one cell's elements (line order) into scr's
  // primitives/conservatives/rhs/density fields. Returns false if the line
  // is unusable (mirrors DiffusionOperator::gatherLine's placeholder guard).
  static bool gatherLine(
    const std::vector<int>& plist,
    ClemParticleContainer::ParticleType* particles,
    ReactionScratch& scr);

  // Fernando-Clem: isobaric closure write-back (rho,e reconstructed at
  // press_held from the reacted T,Y - mirrors DiffusionOperator::scatterLine)
  static void scatterLine(
    ClemParticleContainer::ParticleType* particles,
    const ReactionScratch& scr,
    Diagnostics& diag);
};

} // namespace clem
#endif
