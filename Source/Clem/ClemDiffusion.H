#ifndef CLEMDIFFUSION_H
#define CLEMDIFFUSION_H

#include "PelePhysics.H"

#include "ClemManager.H"
#include "ClemReaction.H"
#include "ClemScratch.H"

namespace clem {

// Fernando-Clem: TransParmType now lives in ClemScratch.H, next to the line
// arena it is always used with - the stirring closure needs it too (it forms
// the line's molecular nu), and it must not have to include this header to
// get it

// Fernando-Clem: molecular diffusion of species and TEMPERATURE along the LEM
// line (module of the ClemAlgorithm composition; see
// InformationImplementation/CouplingLESwithCLEM.md - ClemDiffusion.md is
// deprecated).
//
// Line-internal operator: it moves species mass and energy BETWEEN the elements
// of one cell's line and never between cells (that is the splicing operator).
// The line is PERIODIC (Maxwell's homogeneous sub-grid sample): face n_elem-1
// wraps element n_elem-1 -> 0. Species mass telescopes around the ring, so
// sum_l dm_l dY_k/dt = 0 exactly.
//
// FORMULATION - mass coordinate, CONSTANT PRESSURE. With dm = rho A dx,
// d/dx = A rho d/dm, the species equation
//
//   rho DY_k/Dt = d/dx( Dhat_k dX_k/dx - Y_k sum_j Dhat_j dX_j/dx )
//
// becomes, after the leading rho cancels,
//
//   DY_k/Dt = A^2 d/dm( S_k ),   S_k = F_k - Y_k sum_j F_j,
//   F_k     = (rho Dhat_k) dX_k/dm
//
// and the TEMPERATURE equation (T is the transported thermal variable):
//
//   dT/dt = (A^2/Cp) d/dm( rho lambda dT/dm )     [conduction]
//         + (A^2/Cp) sum_k Cp_k S_k dT/dm         [interspecies transport]
//
// The interspecies term uses the sensible Cp_k, NOT h_k: the formation energy
// is carried by Y and reappears when e = e(rho, T, Y) is reconstructed after
// each T update. Conduction is a flux divergence (telescopes); the interspecies
// term is a pointwise gradient product (does NOT telescope), so the scheme is
// non-conservative in energy - the drift is measured in diffuse()
// (Diagnostics::drift_enthalpy). At constant pressure H = sum(m h), not internal
// energy, is the conserved quantity.
//
// Dhat_k = rho D_k W_k/W is what PelePhysics' transport returns in Ddiag
// (Simple.H: Ddiag[k] = W_k D_k^mix PATM/(RU T) = rho D_k W_k/W), used as the
// coefficient of dX_k/dx; rhoD = rho * Dhat_k is stored. The baro-diffusion term
// (X_k - Y_k) dlnp/dx is dropped: the LEM array is a constant-pressure array.
//
// ISOBARIC CLOSURE - each sub-step holds the element pressure fixed and, after
// updating T and Y, recovers (rho, e) from PYT2RE(p, Y, T): the line breathes
// (rho, hence the element volume mass/rho, changes) while the mass grid is
// fixed. That dilatation (ensemble volume) is what the LES->SGS coupling needs.
//
// DISCRETISATION - central differences on the periodic faces, explicit Euler in
// time. Face coefficients (rho Dhat_k, rho lambda) use a harmonic mean; face
// Y_k and Cp_k are averaged (2nd order).
//
// SUB-STEPPING - the explicit update is conditionally stable while the LES dt is
// not, so the operator advances over N_sub internal sub-steps sized by the
// diffusion-number limit (clem.diffusion_cfl), recomputing the coefficients and
// driving forces from the current (T, Y) at every sub-step. The thermal limit
// uses Cp; the species limit uses rho^2 D_k = rho Dhat_k * W/W_k (the X->Y change
// of variable, see stableSubStep).
class DiffusionOperator
{
public:
  // Fernando-Clem: diagnostics of one diffusion pass (MPI-reduced)
  struct Diagnostics
  {
    // Fernando-Clem: worst sub-step count over the level - the cost knob, and
    // the tell-tale of an under-resolved line (it grows as dm^-2)
    int max_substeps = 0;
    // Fernando-Clem: TOTAL explicit sub-step evaluations over the whole level
    // this LES step = sum over diffused lines of n_sub (MPI ReduceLongSum). This
    // is the real cost of the stage: each unit is one computeRates+advanceLine
    // sweep of a full line. max_substeps only tells the worst line; this tells
    // the work. n_lines is how many lines were diffused, so total/n_lines is the
    // average sub-step count
    amrex::Long total_substeps = 0;
    amrex::Long n_lines = 0;
    // Fernando-Clem: lines that wanted MORE sub-steps than the cap allows and
    // were therefore integrated ABOVE their stability limit. Non-zero means the
    // run is unstable, not merely slow - it is reported loudly for that reason
    amrex::Long n_capped = 0;
    // Fernando-Clem: lines skipped as unusable (empty cell, non-positive mass
    // or density); a healthy run reports 0
    amrex::Long n_skipped = 0;
    // Fernando-Clem: elements where a Y_k came out negative by more than
    // clip_tol and had to be floored + renormalised. Round-off negatives (a
    // species that is identically 0 over half the domain produces them by the
    // bucketful) are floored silently and NOT counted - they would drown the
    // signal. A non-zero count here means a physically significant negative,
    // i.e. the sub-step limit is too loose
    amrex::Long n_clipped = 0;
    // Fernando-Clem: worst |Y_k| that was floored, threshold-free: ~1e-18 is
    // round-off, ~1e-3 is an unstable line
    amrex::Real max_neg_Y = 0.0;

    // Fernando-Clem: conservation drift of the stage (clem.verbose > 1 only -
    // it costs MPI reductions), measured by ClemAlgorithm::advance bracketing
    // conservedTotals() before/after the DiffusionOperator::diffuse() call.
    //
    // DIFFUSION ALONE (clem.do_react=0): the line is closed and the element
    // masses are untouched, so both totals are invariant and both numbers are
    // round-off (~1e-15) unless the scheme is broken.
    //
    // CO-STEPPED WITH REACTION (clem.do_react=1): reaction is interleaved
    // INSIDE diffuse()'s substep loop (see ClemDiffusion.cpp), so this
    // bracket now also spans reaction, and BOTH numbers legitimately stop
    // being round-off - NOT a sign of a broken scheme:
    //   drift_energy: reaction converts chemical (formation) energy to
    //     sensible energy at fixed pressure, so sum(m e) moves by the real
    //     heat release (same reasoning as the pre-existing isobaric-breathing
    //     p dV note below, just a second legitimate contributor now)
    //   drift_species: reaction consumes/produces individual species by
    //     design (that is what combustion means) - only their SUM (total
    //     mass) is conserved, not each mY_k individually, so a large
    //     drift_species under do_react=1 reflects real chemistry, not error.
    //     It remains a meaningful round-off check ONLY when do_react=0.
    //
    // Note the scales. Energy is measured against sum(m |e|), NOT |sum(m e)|:
    // e carries the formation energies, so the signed total is a cancellation
    // and dividing by it inflates a harmless error. Species mass is measured
    // against the TOTAL mass, not against that species' own mass: a radical
    // can hold ~1e-16 g over the whole level, and normalising by that turns
    // clipping noise into a fake 1% error
    amrex::Real drift_energy = 0.0;
    amrex::Real drift_species = 0.0;

    // Fernando-Clem: enthalpy drift of the T-transport, computed inside
    // diffuse() (always, not gated). Constant-pressure diffusion conserves
    // H = sum(m h) on the line, NOT internal energy (the line breathes and does
    // p dV work), so this - not drift_energy - is the meaningful conservation
    // check for the temperature form: |sum(m h)_after - _before| / sum(m |h|).
    // When reaction is co-stepped (react_diag != nullptr), H_after is taken
    // AFTER reactSubstep too - this stays a valid check, not a widened one:
    // the HP-mode reactor is itself adiabatic+isobaric (dH/dt = 0 by
    // construction), so sum(m h) should be invariant across diffusion AND
    // reaction alike, and this number validates BOTH.
    amrex::Real drift_enthalpy = 0.0;

    // Fernando-Clem: stirring (clem.do_stir), counted over the whole level.
    // total_stir_events is the number of triplet maps actually APPLIED;
    // max_stir_events is the worst single line
    amrex::Long total_stir_events = 0;
    int max_stir_events = 0;
    // Fernando-Clem: eddies that were DRAWN (so the clock advanced) but could
    // not be mapped because they spanned fewer than 6 elements. f(l) ~ l^-8/3
    // piles up just above eta, so this is large whenever n_lem < ~6 Delta/eta
    // - and when it is, the realised D_T is far below the target. It is the
    // number that says whether the stage is doing anything
    amrex::Long n_stir_rejected = 0;
    // Fernando-Clem: lines whose stirring sequence was cut short because the
    // substep cap (clem.diffusion_max_substeps) forced the sub-cycle to jump
    // to the end of dt with epochs still pending
    amrex::Long n_stir_truncated = 0;
  };

  // Fernando-Clem: advance every line on the level by dt. tparm is the HOST
  // transport-parameter block (the operator runs on the host, over the particle
  // AoS, like the rest of the Clem module). react_diag, when non-null (i.e.
  // clem.do_diffusion && clem.do_react), makes the substep loop below ALSO
  // call ReactionOperator::reactSubstep after every advanceLine, on the same
  // dt_sub - diffusion and reaction co-stepped over the diffusion-CFL substep
  // (Maxwell [M §3.3.3] steps 14-20), not two independent full-dt stages.
  //
  // nut_strain/filter_width turn STIRRING on (clem.do_stir): the (nu_t, |S|)
  // field built once per step by clem::computeTurbViscosity, and Delta from
  // clem::filterWidth. The sub-cycle then also carries the stirring clock - see
  // the STIRRING block in diffuse(). nullptr/0 leaves the loop exactly as it
  // was, so a non-stirring run is bit-for-bit unchanged
  static Diagnostics diffuse(
    ClemManager& mgr,
    int lev,
    amrex::Real dt,
    const TransParmType* tparm,
    ReactionOperator::Diagnostics* react_diag = nullptr,
    const amrex::MultiFab* nut_strain = nullptr,
    amrex::Real filter_width = 0.0);

private:
  // Fernando-Clem: copy one cell's elements into the arena in line order
  // (sorted by lem_index); false = the line is unusable and must be skipped
  static bool gatherLine(
    LineScratch& lscr,
    const std::vector<int>& plist,
    ClemParticleContainer::ParticleType* particles);

  // Fernando-Clem: write (T, e, Y, rho, vol) back to the particles. The isobaric
  // sub-cycle held press fixed, so it is unchanged; rho breathed, vol = mass/rho
  static void scatterLine(
    const LineScratch& lscr, ClemParticleContainer::ParticleType* particles);

  // Fernando-Clem: (X_k, Cp_k, rho*D_k, rho*lambda, Cp) from the current
  // (rho, T, Y) - the sub-step's coefficients and driving forces
  static void updateProperties(LineScratch& lscr, const TransParmType* tparm);

  // Fernando-Clem: face fluxes (S_k, conduction q, interspecies qint) then the
  // mass-coordinate assembly; fills dYdt/dTdt. Needs updateProperties first
  static void computeRates(LineScratch& lscr, amrex::Real area2);

  // Fernando-Clem: the CFL-stable sub-step from the diffusion-number limit,
  // dt_cfl = cfl * dm_min^2 / max(beta_T, beta_Y) with
  // beta_T = max_l rho lambda A^2 / Cp and beta_Y = max_{l,k} rho^2 D_k A^2
  // (= rho Dhat_k * W/W_k A^2, the X->Y change of variable).
  // Returns a huge value when the line cannot diffuse (no gradients possible)
  static amrex::Real stableSubStep(const LineScratch& lscr, amrex::Real area2);

  // Fernando-Clem: T += dt dTdt, Y += dt dYdt, then (rho, e) from PYT2RE(p, T, Y)
  static void advanceLine(LineScratch& lscr, amrex::Real dt_sub, Diagnostics& diag);

  // Fernando-Clem: THE REARRANGEMENT OPERATOR - one triplet map over the n_seg
  // elements starting at line position `start` (the segment wraps: the line is
  // periodic). ClemStirring says WHEN, HOW BIG and WHERE; this is the map
  // itself.
  //
  // The map lays down three compressed copies of the segment, the middle one
  // reversed, which on an equal-spaced discrete line is exactly the PERMUTATION
  // that takes the elements at offsets 0,3,6,... then 1,4,7,... reversed, then
  // 2,5,8,... For n_seg = 6 that is [0,3,4,1,2,5]. Being a permutation is what
  // makes stirring conservative here: the regrid leaves the line equal-mass, so
  // permuting (rho, T, e, p, Y_k) between elements moves fluid WITHOUT moving
  // mass, and the multiset of element states is preserved exactly. sum(m),
  // sum(m Y_k) and sum(m e) therefore change only by the round-off of
  // re-summing the same values in a different order (~1e-16 relative) - the map
  // itself contributes no conservation error at all, unlike diffusion's
  // non-conservative T-transport. n_seg: positive multiple of 3, <= n_elem
  static void tripletMap(LineScratch& lscr, int start, int n_seg);
};

} // namespace clem
#endif
