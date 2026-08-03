#ifndef CLEMALGORITHM_H
#define CLEMALGORITHM_H

#include "ClemDiffusion.H"
#include "ClemManager.H"
#include "ClemReaction.H"
#include "ClemSplicing.H"
#include "ClemStirring.H"

namespace clem {

// Fernando-Clem: per-step subgrid algorithm (alpaca-style module): the stage
// sequence of one CLEM step is assembled here from the operator modules, not
// inside PeleC. Current composition (Maxwell chain):
//   1. scale the accumulated face mass fluxes (time centering)
//   2. LesCoupling::isentropicPressureUpdate - LES->SGS pressure work
//   3. Regridder - restore n_lem equal-mass elements per cell
//   4-5. DiffusionOperator + ReactionOperator (clem.do_diffusion,
//      clem.do_react), CO-STEPPED when both are on: Maxwell [M §3.3.3] steps
//      14-20 - Delta t_diff is computed from the diffusion CFL, diffusion and
//      reaction are integrated TOGETHER over that substep (reaction right
//      after each diffusion substep, external source forced to zero - see
//      ClemReaction.H), repeated until the full LES dt is covered. This is
//      SEQUENTIAL (Lie) splitting at the SUBSTEP granularity, not two
//      independent full-dt stages. If diffusion is off, Delta t_diff
//      degenerates to the full LES dt and reaction alone runs once
//      (ReactionOperator::react).
//   6. STIRRING (clem.do_stir) - NOT a stage of its own. A triplet map is an
//      instantaneous rearrangement, so it rides INSIDE the 4-5 sub-cycle as a
//      third clip on Delta t_diff: the (nu_t, |S|) field is built here, once
//      per step, from the resolved LES state, and the per-line clock lives in
//      DiffusionOperator::diffuse. See the STIRRING block there.
//   7. SplicingOperator - inter-cell transport with the LES mass fluxes
//   8. LesCoupling::writeBackToLes - SGS->LES filtered (rho Y_k, rho e)
//
// Diffusion sits AFTER the regrid on purpose: the regrid is what restores the
// equal-mass line the mass-coordinate discretisation is written for, and it is
// also what merges slivers - running it after diffusion would smear the
// gradients diffusion just resolved.
class Algorithm
{
public:
  // Fernando-Clem: diagnostics of one subgrid step
  struct StepDiagnostics
  {
    SpliceBalance balance;
    amrex::Real total_mass = 0.0;
    // Fernando-Clem: max|rho_filtered - rho_LES|/rho_LES (coupling quality)
    amrex::Real coupling_rho_err = 0.0;
    DiffusionOperator::Diagnostics diffusion;
    ReactionOperator::Diagnostics reaction;
  };

  // Fernando-Clem: advance the subgrid by dt. flux_scale converts the
  // accumulated MOL-stage fluxes to the time-centered flux (0.5 for the
  // two-stage Heun scheme); state_bc must be FillPatched with >= 1 ghost
  // cell (bcnormal inflow state); s_new is the LES state at t^{n+1} (its
  // pressure drives the isentropic stage; its composition/energy are
  // overwritten by the write-back when clem.couple_back is on); nstep drives
  // the periodic particle dumps; tparm is PeleC's HOST transport-parameter
  // block, needed by the diffusion stage; les is a read-only view of the LES
  // model PeleC is actually running, and is what makes the stirring rate a
  // property of the resolved field rather than a prescribed number
  static StepDiagnostics advance(
    ClemManager& mgr,
    int lev,
    amrex::Real dt,
    const amrex::MultiFab& state_bc,
    amrex::MultiFab& s_new,
    amrex::MultiFab& s_old,
    int rho_indx,
    int mom_indx,
    int eden_indx,
    int temp_indx,
    int eint_indx,
    int spec_indx,
    amrex::Real flux_scale,
    int nstep,
    const TransParmType* tparm,
    const LesView& les = LesView{});
};

} // namespace clem
#endif
