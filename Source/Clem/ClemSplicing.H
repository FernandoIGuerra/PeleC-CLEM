#ifndef CLEMSPLICING_H
#define CLEMSPLICING_H

#include "ClemManager.H"

namespace clem {

// Fernando-Clem: global mass-balance record of one splicing pass
// (MPI-reduced); dM_total must equal mass_in - mass_out
struct SpliceBalance
{
  amrex::Real mass_in = 0.0;  // injected through inflow boundaries
  amrex::Real mass_out = 0.0; // removed through outflow boundaries
};

// Fernando-Clem: splicing/advection operator (Maxwell thesis, steps 23-24).
// Moves LEM elements between LES cells according to the signed face mass
// fluxes in ClemManager::fluxFaces() (positive = incoming, negative =
// outgoing, [g/s]; multiplied by dt here). The arrival ordering strategy is
// the compile-time config::splice_policy (side_aware vs append_downstream)
// so both can be compared. Requires state_bc FillPatched with >= 1 ghost
// cell: ghost values supply the incoming-fluid state at inflow boundaries.
class SplicingOperator
{
public:
  static SpliceBalance splice(
    ClemManager& mgr,
    int lev,
    amrex::Real dt,
    const amrex::MultiFab& state_bc,
    int rho_indx,
    int temp_indx,
    int eint_indx,
    int spec_indx);

private:
  // Fernando-Clem: per-cell face walk - full moves, splits, boundary
  // injection/removal; returns the local (unreduced) balance
  static SpliceBalance advectCells(
    ClemManager& mgr,
    int lev,
    amrex::Real dt,
    const amrex::MultiFab& state_bc,
    int rho_indx,
    int temp_indx,
    int eint_indx,
    int spec_indx);

  // Fernando-Clem: Maxwell step 24 - re-index each line so arrivals sit on
  // the end matching their entry face (left entries -> lowest indices,
  // residents keep order, right/transverse -> highest); clears the transient
  // entry_side/xflux bookkeeping
  static void reindexAfterSplicing(ClemManager& mgr, int lev);
};

} // namespace clem
#endif
