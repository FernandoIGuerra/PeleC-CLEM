#include <cmath>

#include <AMReX_Utility.H>

#include "ClemAlgorithm.H"
#include "ClemLesCoupling.H"
#include "ClemReaction.H"
#include "ClemRegrid.H"

namespace clem {

namespace {

// Fernando-Clem: the totals the diffusion stage must leave alone. It only
// redistributes energy and species mass WITHIN a line (closed ends, element
// masses untouched), so every one of these is invariant across the stage
struct ConservedTotals
{
  amrex::Real e = 0.0;               // sum(m e) - total internal energy
  amrex::Real abs_e = 0.0;           // sum(m |e|) - the SCALE of that total
  amrex::Real mass = 0.0;            // sum(m)
  amrex::Real mY[NUM_SPECIES] = {0.0}; // sum(m Y_k) - mass of each species
};

// Fernando-Clem: costs 2 + NUM_SPECIES MPI reductions, hence the verbose gate
ConservedTotals
conservedTotals(const ClemManager& mgr, const int lev)
{
  ConservedTotals t;
  t.e           = mgr.sumMassWeighted(lev, RealData::eint);
  t.abs_e       = mgr.sumMassWeightedAbs(lev, RealData::eint);
  t.mass        = mgr.sumParticleReal(lev, RealData::mass);
  
  for (int k = 0; k < NUM_SPECIES; ++k) {
    t.mY[k]     = mgr.sumMassWeighted(lev, RealData::Y0 + k);
  }
  return t;
}

// Fernando-Clem: energy drift against sum(m |e|), NOT against |sum(m e)|: e
// carries the species formation energies, so the signed total is a sum of
// opposite-sign contributions and can sit near zero by cancellation - dividing
// by it reports a fake error
amrex::Real
driftEnergy(const ConservedTotals& t0, const ConservedTotals& t1)
{
  return (t0.abs_e > 0.0) ? std::abs(t1.e - t0.e) / t0.abs_e : 0.0;
}

// Fernando-Clem: species drift against the TOTAL mass, NOT against each
// species' own mass: a radical holds ~1e-16 g over the whole level, and
// normalising by that turns round-off into a fake percent-level error. Against
// the total mass the number means what it should - the fraction of the mass in
// the system that the stage misplaced
amrex::Real
driftSpecies(const ConservedTotals& t0, const ConservedTotals& t1)
{
  if (t0.mass <= 0.0) {
    return 0.0;
  }
  amrex::Real drift = 0.0;
  for (int k = 0; k < NUM_SPECIES; ++k) {
    drift = amrex::max<amrex::Real>(
      drift, std::abs(t1.mY[k] - t0.mY[k]) / t0.mass);
  }
  return drift;
}

} // namespace

Algorithm::StepDiagnostics
Algorithm::advance(
  ClemManager& mgr,
  const int lev,
  const amrex::Real dt,
  const amrex::MultiFab& state_bc,
  amrex::MultiFab& s_new,
  amrex::MultiFab& s_old,
  const int rho_indx,
  const int mom_indx,
  const int eden_indx,
  const int temp_indx,
  const int eint_indx,
  const int spec_indx,
  const amrex::Real flux_scale,
  const int nstep,
  const TransParmType* tparm)
{
  BL_PROFILE("clem::Algorithm::advance()");
  AMREX_ALWAYS_ASSERT(mgr.isDefined());

  // Fernando-Clem: stage - time centering of the accumulated face fluxes
  if (flux_scale != 1.0) {
    mgr.fluxFaces().mult(flux_scale);
  }

  // Fernando-Clem: write-back mode. freeze_species=1 -> Maxwell conserved-delta +
  // projection (subgrid owns species; the LES did NOT advect rho Y_k); needs the
  // CONSERVED filtered species density at the START of the step (before any
  // subgrid stage) and after splicing. freeze_species=0 -> React.cpp increment.
  // The SAME conserved pre/post snapshots feed the LES-vs-SGS consistency
  // instrument (clem.verbose > 0): per-step d(rhoY_k)/d(rhoE) of the CLEM
  // machinery against the LES change, split by channel.
  const bool couple_back = ClemManager::coupleBack();
  const bool maxwell = ClemManager::freezeSpecies();
  const bool do_diag = (ClemManager::verbose() > 0);
  amrex::MultiFab cons_pre, cons_post;
  if (do_diag || (couple_back && maxwell)) {
    cons_pre.define(s_new.boxArray(), s_new.DistributionMap(), NUM_SPECIES + 1, 0);
    LesCoupling::snapshotFilteredMean(mgr, lev, cons_pre, /*conserved=*/true);
  }

  // Fernando-Clem: stage - LES -> SGS pressure coupling: every element is
  // driven from its own pressure p_k to the resolved cell pressure P(s_new)
  // as isentropic work (clem.pressure_coupling)
  //Comment so far the pressure is updates using the isentropic relations in temperature and volumen.
  //Also the energy accounts for this change
  if (ClemManager::pressureCoupling()) {
    LesCoupling::isentropicPressureUpdate( mgr, lev, s_new, rho_indx, temp_indx, spec_indx);
  }

  // Fernando-Clem: stage - regrid back to n_lem equal-mass elements
  //regridding must be every step MANDATORY
  Regridder::regrid(mgr, lev);


  StepDiagnostics diag;

  // Fernando-Clem: coupling increment brackets. sgs_incr first holds the PRE
  // filtered mean (taken here, AFTER the isentropic update and regrid so their
  // bookkeeping is excluded) and BEFORE the process stages; after diffusion/
  // reaction and BEFORE splice we subtract it from the POST mean so it holds
  // the process-only increment (Ybar,ebar). Transport (splice) never enters it.
  amrex::MultiFab sgs_incr;
  amrex::MultiFab sgs_post;
  if (couple_back && !maxwell) {
    sgs_incr.define(s_new.boxArray(), s_new.DistributionMap(), NUM_SPECIES + 1, 0);
    sgs_post.define(s_new.boxArray(), s_new.DistributionMap(), NUM_SPECIES + 1, 0);
    LesCoupling::snapshotFilteredMean(mgr, lev, sgs_incr); // pre
  }

  // Fernando-Clem: stage - molecular diffusion along the regridded line
  // (species + temperature). Line-internal: it redistributes species mass and
  // energy among the elements of a cell and leaves the element masses alone.
  // The Delta t_diff sub-cycle itself (t_clem = 0; while (t_clem < dt):
  // recompute Delta t_diff from the CFL, clip to what remains, diffuse+react
  // together over it, t_clem += Delta t_diff) is implemented per-line inside
  // DiffusionOperator::diffuse (ClemDiffusion.cpp) - see the comment just
  // below for why it is composed here rather than unrolled at this level.
  if (ClemManager::doDiffusion()) {
    const bool check_conservation = (ClemManager::verbose() > 1);
    ConservedTotals before;
    if (check_conservation) {
      before = conservedTotals(mgr, lev);
    }

    if (ClemManager::doReact()) {
      // Fernando-Clem: stage - diffusion and reaction CO-STEPPED over the
      // diffusion-CFL substep (Maxwell [M §3.3.3] steps 14-20: compute
      // Delta t_diff from the CFL restriction, diffuse+react TOGETHER over
      // that substep, advance, repeat until dt is covered). NOT two
      // independent full-dt stages (diffusion sub-cycling internally for the
      // whole dt, then reaction running once on the fully-diffused end
      // state) - Delta t_diff bounds both operators together, and reaction
      // runs right after EACH diffusion substep, inside diffuse()'s own loop
      // (see ReactionOperator::reactSubstep, called from ClemDiffusion.cpp).
      diag.diffusion = DiffusionOperator::diffuse(mgr, lev, dt, tparm, &diag.reaction);
    } else {
      diag.diffusion = DiffusionOperator::diffuse(mgr, lev, dt, tparm);
    }

    if (check_conservation) {
      const ConservedTotals after     = conservedTotals(mgr, lev);
      diag.diffusion.drift_energy     = driftEnergy(before, after);
      diag.diffusion.drift_species    = driftSpecies(before, after);
    }
  } else if (ClemManager::doReact()) {
    // Fernando-Clem: stage - subgrid chemistry alone. With diffusion off
    // there is no diffusion CFL to bound a substep, so Maxwell's
    // Delta t_diff degenerates to the full LES dt: react once, on the whole
    // step, directly against the particles (ReactionOperator::react).
    // SEQUENTIAL (Lie) splitting either way - the reactor's external source
    // is always forced to zero (see ClemReaction.H for why).
    diag.reaction = ReactionOperator::react(mgr, lev, dt);
  }

  // Fernando-Clem: (future stage) stirring is composed here

  // Fernando-Clem: close the increment bracket BEFORE splice - sgs_incr becomes
  // (POST filtered mean) - (PRE), i.e. exactly what diffusion+reaction changed
  // the filtered composition/energy by, at the same (pre-splice) cell config
  
  //if (couple_back) {
  //  LesCoupling::snapshotFilteredMean(mgr, lev, sgs_post); // post
  //  amrex::MultiFab::Subtract(sgs_incr, sgs_post, 0, 0, NUM_SPECIES + 1, 0);
  //  sgs_incr.mult(-1.0, 0, NUM_SPECIES + 1, 0); // incr = post - pre
  //}

  // Fernando-Clem: stage - splicing (inter-cell transport)
  diag.balance = SplicingOperator::splice(mgr, lev, dt, state_bc, rho_indx, temp_indx, eint_indx, spec_indx);

  // Fernando-Clem: ensemble-volume renormalization (clem.vol_renorm) - put the
  // line back into the fixed-volume cell AFTER splicing and BEFORE it is
  // filtered/written back. Catches ALL upstream volume defects (isentropic
  // work, isobaric diffusion breathing, regrid + splice quantization). The p dV
  // of the line's breathing is NOT seen by the fixed-volume LES cell, and it
  // polluted the filtered (rho e) by ~ rho*e * vol_defect (8x the LES per-step
  // energy change with diffusion on). Isentropic uniform rescale to
  // sum(vol_p) = V_cell converts the breathing work back into internal energy
  // (net ~ constant-volume closure of the isobaric stages).
  amrex::Real vol_renorm_corr = 0.0;
  if (ClemManager::volRenorm()) {
    vol_renorm_corr = LesCoupling::renormalizeEnsembleVolume(mgr, lev);
  }

  // Fernando-Clem: close the conserved bracket AFTER splicing, BEFORE the
  // write-back (so s_new still holds the pure LES change). Shared by the
  // LES-vs-SGS instrument and the Maxwell write-back below.
  if (do_diag || (couple_back && maxwell)) {
    cons_post.define(s_new.boxArray(), s_new.DistributionMap(), NUM_SPECIES + 1, 0);
    LesCoupling::snapshotFilteredMean(mgr, lev, cons_post, /*conserved=*/true);
  }

  // Fernando-Clem: LES-vs-SGS conserved instrument - per-step effect of the two
  // machineries on (rho Y_k) and (rho e), split by channel + accumulated sync
  // offset. Run in the LES-driven mode (freeze=0, LES = truth) to VALIDATE the
  // CLEM energy and species-mass bookkeeping channel by channel.
  if (do_diag) {
    LesCoupling::conservedConsistency(
      mgr, lev, cons_pre, cons_post, s_new, s_old, rho_indx, eint_indx,
      spec_indx);
  }

  // Fernando-Clem: stage - SGS -> LES species write-back.
  if (couple_back && maxwell) {
    // Maxwell conserved-delta + projection: feed back (rho Y_k)^n_LES + the
    // CLEM conserved delta and project onto sum_k rho Y_k = rho_LES (see
    // writeBackProjected).
    diag.coupling_rho_err = LesCoupling::writeBackProjected(
                            mgr, lev, s_new, s_old, cons_pre, cons_post,
                            rho_indx, mom_indx, eden_indx, temp_indx, eint_indx,
                            spec_indx);
  } else if (couple_back) {
    // React.cpp-style INCREMENT feedback (freeze_species=d0): the LES keeps its
    // own resolved transport; only the subgrid process increment is added.
    LesCoupling::snapshotFilteredMean(mgr, lev, sgs_post); // post
    amrex::MultiFab::Subtract(sgs_incr, sgs_post, 0, 0, NUM_SPECIES + 1, 0);
    sgs_incr.mult(-1.0, 0, NUM_SPECIES + 1, 0); // incr = post - pre
    diag.coupling_rho_err = LesCoupling::writeBackToLes(
                                              mgr, lev, s_new, s_old, sgs_incr, rho_indx, mom_indx, eden_indx, temp_indx,
                                              eint_indx, spec_indx);
  }

  //###########################################################################################################
  //                              DIAGNOSTICS
  //###########################################################################################################
  diag.total_mass = mgr.sumParticleReal(lev, RealData::mass);
  if (ClemManager::verbose() > 0) {
    // Fernando-Clem: volume defect = max|sum(vol_p)/V_cell - 1|; tests the
    // consistency of the ensemble volume with the cell (density hypothesis)
    const amrex::Real vol_defect = mgr.maxCellVolumeDefect(lev);
    amrex::Print() << "CLEM advance: mass = " << diag.total_mass
                   << ", in = " << diag.balance.mass_in
                   << ", out = " << diag.balance.mass_out
                   << ", rho coupling err = " << diag.coupling_rho_err
                   << ", vol defect = " << vol_defect
                   << ", vol renorm corr = " << vol_renorm_corr << '\n';
    if (ClemManager::doDiffusion()) {
      const amrex::Long nl = diag.diffusion.n_lines;
      const amrex::Real avg_sub =
        (nl > 0) ? static_cast<amrex::Real>(diag.diffusion.total_substeps) /
                     static_cast<amrex::Real>(nl)
                 : 0.0;
      amrex::Print() << "CLEM diffusion: substep calls = "
                     << diag.diffusion.total_substeps << " over " << nl
                     << " lines (avg = " << avg_sub
                     << ", max = " << diag.diffusion.max_substeps << ")"
                     << " (capped lines = " << diag.diffusion.n_capped << ")"
                     << ", skipped lines = " << diag.diffusion.n_skipped
                     << ", clipped elements = " << diag.diffusion.n_clipped
                     << ", max neg Y = " << diag.diffusion.max_neg_Y
                     << ", drift(h) = " << diag.diffusion.drift_enthalpy;
      if (ClemManager::verbose() > 1) {
        // Fernando-Clem: drift(e) reflects the PHYSICAL p dV work of the
        // isobaric line (rho breathes), not a numerical error - drift(h)
        // above is the conservation check for the constant-pressure
        // T-transport. With clem.do_react=1, reaction is co-stepped INSIDE
        // diffuse() (see ClemDiffusion.H's Diagnostics comment), so both
        // drift(e) and drift(Y) also carry reaction's real heat release /
        // species conversion here - large numbers are expected, not a sign
        // of a broken scheme. drift(Y) is a round-off-only sanity check
        // ONLY when clem.do_react=0.
        amrex::Print() << ", drift(e) = " << diag.diffusion.drift_energy
                       << ", drift(Y) = " << diag.diffusion.drift_species;
      }
      amrex::Print() << '\n';
    }
    if (ClemManager::doReact()) {
      amrex::Print() << "CLEM reaction: lines = " << diag.reaction.n_lines
                     << ", skipped lines = " << diag.reaction.n_skipped
                     << ", max neg Y = " << diag.reaction.max_neg_Y
                     << ", retried lines = " << diag.reaction.n_retried_lines
                     << " (max retries = " << diag.reaction.max_retries << ")"
                     << '\n';
    }
  }

  // Fernando-Clem: periodic particle dumps for the transport validation
  const int plot_int = ClemManager::plotInt();
  if (plot_int > 0 && nstep % plot_int == 0) {
    mgr.writePlotFile(amrex::Concatenate("plt_clem_", nstep, 5), "particles");
  }

  return diag;
}

} // namespace clem
