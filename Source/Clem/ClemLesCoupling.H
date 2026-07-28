#ifndef CLEMLESCOUPLING_H
#define CLEMLESCOUPLING_H

#include "ClemManager.H"

namespace clem {

// Fernando-Clem: the two coupling operators between the LES (super-grid) and
// the LEM subgrid (module of the ClemAlgorithm composition):
//
//   LES -> SGS: isentropicPressureUpdate - every element is driven from ITS
//   OWN stored pressure p_k (Maxwell element state) to the new resolved cell
//   pressure P(S_new) as isentropic work (composition frozen):
//     T/T_old = (P/p_k)^((g-1)/g),  rho/rho_old = (P/p_k)^(1/g)
//   with the local gamma = cp/cv, and the energy updated in DIFFERENCE form
//   e += e(T_new,Y) - e(T_old,Y) so a constant T-representation bias cannot
//   ratchet energy. (Historical note kept from the old code: if a per-step
//   energy pump reappears at the front, the suspect is the stored-pressure
//   reference interacting with the write-back correction; the alternative
//   reference is the cell P(S_old).)
//
//   SGS -> LES: writeBackToLes - React.cpp-style INCREMENT feedback. Instead of
//   overwriting the LES composition with the filtered line (which replaces the
//   2nd-order LES transport with the 1st-order spliced transport and, in a
//   multi-species EOS, shifts pressure and puts a spurious velocity bulge at
//   fronts), the LES keeps its own resolved transport and the subgrid feeds
//   back only the change its PROCESSES (diffusion + reaction) made to the
//   filtered composition/energy over the step:
//     rhoY_k += rho_LES * ( Ybar_k^post - Ybar_k^pre )
//   The increment is bracketed between "after regrid" (pre) and "before splice"
//   (post), so transport (splice) and the isentropic/regrid bookkeeping never
//   reach the LES. In pure advection there are no processes -> increment is
//   exactly zero -> the LES stays bit-identical to stock PeleC (no bulge).
//   sum_k of the increment is zero (mass conserved by the processes), so
//   sum_k rhoY_k stays = rho_LES. clem.couple_back_energy=1 also feeds back the
//   filtered internal-energy increment (rho e += rho_LES * (ebar^post-ebar^pre)).
class LesCoupling
{
public:
  static void isentropicPressureUpdate(
    ClemManager& mgr,
    int lev,
    const amrex::MultiFab& s_new,
    int rho_indx,
    int temp_indx,
    int spec_indx);

  // Fernando-Clem: ENSEMBLE-VOLUME renormalization (clem.vol_renorm). The
  // isobaric line closure lets the line breathe (sum vol_p drifts from V_cell);
  // the p dV work of that breathing pollutes the filtered (rho e) because the
  // LES cell has FIXED volume (measured: energy gap ~ rho*e * vol_defect, 8x
  // the LES per-step energy change with diffusion on). This puts the gas back
  // in the box: per cell, every element is compressed/expanded ISENTROPICALLY
  // by the uniform ratio f = V_cell / sum(vol_p), converting the breathing work
  // back into internal energy (net: isobaric stage + renorm ~ constant-volume
  // process). Multi-species analog of Maxwell's perfect-gas identity
  // "equal p <=> equal rho e" which does NOT hold for a real gas.
  // Returns the MPI-reduced max|f - 1| (the defect just corrected).
  static amrex::Real renormalizeEnsembleVolume(ClemManager& mgr, int lev);

  // Fernando-Clem: filtered per-cell mean of the current subgrid line.
  //   conserved=false (default): out(iv,k) = Ybar_k = sum(m Y_k)/sum(m),
  //     out(iv,NUM_SPECIES) = ebar = sum(m e)/sum(m)   (Favre RATIO form)
  //   conserved=true: out(iv,k) = sum(m Y_k)/V_cell = (rho Y_k)_filtered,
  //     out(iv,NUM_SPECIES) = sum(m e)/V_cell = (rho e)_filtered  (CONSERVED)
  // Empty cells stay 0. out must have NUM_SPECIES+1 components.
  static void snapshotFilteredMean(
    ClemManager& mgr,
    int lev,
    amrex::MultiFab& out,
    bool conserved = false);

  // Fernando-Clem: apply the pre-computed filtered increment sgs_incr
  // (NUM_SPECIES+1 comps, = post - pre from snapshotFilteredMean) to the LES
  // state. Returns the coupling-quality metric max|rho_filtered - rho_LES|/
  // rho_LES over the level (MPI-reduced) - zero up to splicing round-off when
  // the subgrid transports exactly the LES mass fluxes
  static amrex::Real writeBackToLes(
    ClemManager& mgr,
    int lev,
    amrex::MultiFab& s_new,
    amrex::MultiFab& s_old,
    const amrex::MultiFab& sgs_incr,
    int rho_indx,
    int mom_indx,
    int eden_indx,
    int temp_indx,
    int eint_indx,
    int spec_indx);

  // Fernando-Clem: LES-vs-SGS CONSERVED consistency instrument. Per cell, per
  // step, compares the EFFECT of the two machineries on the same conserved
  // filtered quantities:
  //   species: d(rho Y_k)_CLEM = cons_post_k - cons_pre_k   (splice + diffusion
  //            + reaction on the line, conserved form sum(m Y)/V_cell)
  //     vs     d(rho Y_k)_LES  = s_new - s_old (rho Y_k)     (Godunov flux
  //            [+ LES diffusion])
  //   energy:  d(rho e)_CLEM   = cons_post[NS] - cons_pre[NS] (isentropic work
  //            + spliced transport [+ line diffusion])
  //     vs     d(rho e)_LES    = s_new - s_old (rho e from eint_indx)
  // Also reports the accumulated SYNC OFFSET |cons_pre - s_old| for both
  // channels (how far the subgrid has drifted from the LES baseline).
  // Run with the LES-driven mode (freeze_species=0, couple_back=0/increment) so
  // the LES is the trusted truth; the gaps then localize the CLEM bookkeeping
  // error by channel. cons_pre/post from snapshotFilteredMean(conserved=true),
  // pre = start of the CLEM step, post = after splicing, both before write-back.
  // Returns the MPI-reduced max species gap.
  static amrex::Real conservedConsistency(
    ClemManager& mgr,
    int lev,
    const amrex::MultiFab& cons_pre,
    const amrex::MultiFab& cons_post,
    const amrex::MultiFab& s_new,
    const amrex::MultiFab& s_old,
    int rho_indx,
    int eint_indx,
    int spec_indx);

  // Fernando-Clem: Maxwell-faithful CONSERVED-DELTA + PROJECTION species write-
  // back (freeze_species=1 path). Feeds back the subgrid's conserved species
  // change over the CLEM step and enforces mass realizability against the LES
  // density:
  //   (rho Y_k)^{n+1} = (rho Y_k)^n_LES + [ cons_post_k - cons_pre_k ]      (add
  //                     the CLEM conserved delta, incl. splice transport)
  //   then PROJECT so sum_k (rho Y_k)^{n+1} = rho_LES^{n+1} exactly.
  // cons_pre/cons_post are conserved-form snapshots (snapshotFilteredMean with
  // conserved=true) at the start of the step and after splicing.
  // clem.couple_back_energy selects the energy treatment (see ClemManager.H):
  //   1 = conserved delta: (rho e)^{n+1} = (rho e)^n_LES + [cons_post -
  //       cons_pre]_(rho e) (CLEM isentropic work + spliced transport).
  //   2 = keep-T reconciliation: rho e := rho * e(T_pre, Y_new) with T_pre =
  //       the pre-write-back LES temperature (temp_indx, computeTemp-consistent)
  //       so the composition write-back cannot shift T - the conduction
  //       operator then never sees the EOS T-kick of dY (the velocity-deficit
  //       mechanism with diffusion on). Injects the formation-energy of dY at
  //       fixed T (non-conservative; monitor for pumping).
  // rho E is rebuilt with the LES kinetic energy. Returns the MPI-reduced
  // max|projection factor - 1| (should be ~rho_err).
  static amrex::Real writeBackProjected(
    ClemManager& mgr,
    int lev,
    amrex::MultiFab& s_new,
    const amrex::MultiFab& s_old,
    const amrex::MultiFab& cons_pre,
    const amrex::MultiFab& cons_post,
    int rho_indx,
    int mom_indx,
    int eden_indx,
    int temp_indx,
    int eint_indx,
    int spec_indx);
};

} // namespace clem
#endif
