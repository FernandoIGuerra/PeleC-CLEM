#ifndef CLEMMANAGER_H
#define CLEMMANAGER_H

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <AMReX_BaseFab.H>
#include <AMReX_BCRec.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include "ClemParticleContainer.H"

namespace clem {

// Fernando-Clem: ClemManager owns the particle container and the cell->particle
// map; PeleC talks only to this class (unlike spray, which smears container +
// lifecycle across PeleC statics)
class ClemManager
{
public:
  // Fernando-Clem: per-grid fab of particle-index lists (host side). The
  // interface is fixed so a GPU-friendly CSR/bin layout can replace this
  // storage later without touching callers
  using CellVectors = std::map<int, amrex::BaseFab<std::vector<int>>>;

  ClemManager() = default;

  // Fernando-Clem: runtime params from "clem.*" (verbosity, freeze flag,
  // plot interval; layout and n_lem are compile time, see ClemIndex.H)
  static void readParams();

  // Fernando-Clem: freeze rho*Y_k in the LES advance so the subgrid owns the
  // species transport (clem.freeze_species, default off - study switch)
  static bool freezeSpecies() { return m_freeze_species; }

  // Fernando-Clem: particle plotfile interval in steps (clem.plot_int,
  // <= 0 disables)
  static int plotInt() { return m_plot_int; }
  static int verbose() { return m_verbose; }

  // Fernando-Clem: coupling stage switches for the strategy study
  // (clem.pressure_coupling: LES->SGS isentropic update;
  //  clem.couple_back: SGS->LES filtered composition write-back;
  //  clem.couple_back_energy: ALSO replace the LES internal energy with the
  //  filtered particle energy - 1st-order spliced transport overwrites the
  //  2nd-order LES energy, known to produce a velocity dipole at fronts;
  //  kept only for comparison. Energy should enter via the reaction source)
  static bool pressureCoupling() { return m_pressure_coupling; }
  static bool coupleBack() { return m_couple_back; }
  static bool coupleBackEnergy() { return m_couple_back_energy != 0; }
  // Fernando-Clem: energy write-back mode (clem.couple_back_energy):
  //   0 = off (energy untouched; computeTemp shifts T when Y changes)
  //   1 = conserved delta: rho e += [cons_post - cons_re]_(rho e) (measured
  //       neutral in advection)
  //   2 = keep-T reconciliation: rho e := rho * e(T_pre, Y_new) - the LES
  //       temperature is PRESERVED under the composition write-back, so the
  //       conduction operator never sees the EOS T-shift of dY (the -30 cm/s
  //       deficit mechanism). Non-conservative by the formation-energy of dY
  //       per step - watch the instrument for pumping.
  static int coupleBackEnergyMode() { return m_couple_back_energy; }

  // Fernando-Clem: ensemble-volume renormalization (clem.vol_renorm, default
  // on). After the subgrid process stages, each cell's line is isentropically
  // rescaled so sum(vol_p) = V_cell exactly. Removes the p dV bookkeeping of
  // the ISOBARIC line closure (the line "breathes" while the LES cell has fixed
  // volume) - measured as the energy-channel gap ~ rho*e * vol_defect that
  // collapses the coupling when diffusion is on
  static bool volRenorm() { return m_vol_renorm; }

  // Fernando-Clem: regrid cadence in steps (clem.regrid_int, default 1 =
  // every step). Less frequent regridding reduces the sliver-merge
  // dilution (artificial Y diffusion at the element scale)
  static int regridInt() { return m_regrid_int; }

  // Fernando-Clem: molecular diffusion along the LEM line
  // (clem.do_diffusion). When the subgrid owns diffusion, the LES must not
  // also diffuse the same species/energy - see the note in ClemDiffusion.H
  static bool doDiffusion() { return m_do_diffusion; }

  // Fernando-Clem: diffusion-number limit of the explicit line update
  // (clem.diffusion_cfl, default 1/2 = the 2nd-order-stencil limit); it sets
  // the internal sub-step count, so lowering it costs time linearly
  static amrex::Real diffusionCfl() { return m_diffusion_cfl; }

  // Fernando-Clem: hard cap on the internal sub-steps of one line
  // (clem.diffusion_max_substeps). A runaway count means the line is
  // under-resolved (the limit scales as dm^-2), not that the cap is too low;
  // the cap only keeps a pathological cell from stalling the run
  static int diffusionMaxSubsteps() { return m_diffusion_max_substeps; }

  // Fernando-Clem: subgrid chemistry along the LEM line (clem.do_react).
  // SEQUENTIAL (Lie) operator splitting: reaction runs AFTER diffusion, over
  // the full dt, with the reactor's external source set to zero - diffusion
  // has already been applied to the particles, so reaction just integrates
  // chemistry on the post-diffusion state. Chosen over coupling diffusion's
  // rate into CVODE's rhs because the particle state is PRIMITIVE (T, Y_k),
  // and building a consistent CONSERVATIVE rhs from a primitive line would be
  // its own source of error for no accuracy benefit at dt this small (see
  // ClemReaction.H). When the subgrid owns reaction, the LES must not also
  // react the same species/energy (mirrors clem.do_diffusion's contract)
  static bool doReact() { return m_do_react; }

  // Fernando-Clem: max over cells of |sum(vol_p)/V_cell - 1| (MPI-reduced) -
  // the ensemble volume defect; grows if pressure work and splicing leave
  // the line volume inconsistent with the cell
  amrex::Real maxCellVolumeDefect(int lev) const;

  // Fernando-Clem: create the container from the AMR hierarchy (production
  // path, mirrors SprayPC creation)
  void define(amrex::AmrCore* amr);

  // Fernando-Clem: single-level define for unit tests
  void defineSingleLevel(
    const amrex::Geometry& geom,
    const amrex::DistributionMapping& dmap,
    const amrex::BoxArray& ba);

  // Fernando-Clem: place config::n_lem elements per cell, evenly spaced along
  // x; state values are placeholders until the coupling fills them from the
  // LES state
  void initParticles(int lev);

  // Fernando-Clem: t=0 fill - every element takes the state of its host LES
  // cell; vol = V_cell/n_lem and mass = rho*vol so summing the elements
  // recovers the cell (consistent with the filter definition)
  void setParticlesFromState(
    int lev,
    const amrex::MultiFab& state,
    int rho_indx,
    int temp_indx,
    int eint_indx,
    int spec_indx);

  // Fernando-Clem: (re)build the cell->particle map; must be called after any
  // Redistribute/regrid before the map is used
  void rebuildCellMap(int lev);
  void invalidateCellMap() { m_cell_map_valid = false; }
  bool IsCellMapValid() const { return m_cell_map_valid; }

  // Fernando-Clem: particle indices (into the grid AoS) living in cell iv
  const std::vector<int>& particlesInCell(int lev, int grid_id, const amrex::IntVect& iv) const;

  // Fernando-Clem: structural check - every cell owns exactly config::n_lem
  // elements and every particle is mapped exactly once (reduced over ranks)
  bool verifyCellMap(int lev) const;

  // Fernando-Clem: wraps Redistribute and invalidates the map so a stale map
  // cannot be used silently
  void redistribute(int lbase = 0);

  void writePlotFile(const std::string& dir, const std::string& name) const;

  ClemParticleContainer& particleContainer() { return *m_pc; }
  const ClemParticleContainer& particleContainer() const { return *m_pc; }
  bool isDefined() const { return m_pc != nullptr; }

  amrex::Long totalParticleCount(int lev, bool only_local = false) const;

  // Fernando-Clem: MPI-reduced sum of one particle real component (used for
  // conservation checks against the LES state)
  amrex::Real sumParticleReal(int lev, int comp) const;

  // Fernando-Clem: MPI-reduced sum of mass * rdata(comp) - conserved totals
  // of specific quantities (e.g. comp = eint gives total internal energy)
  amrex::Real sumMassWeighted(int lev, int comp) const;

  // Fernando-Clem: same, but of |rdata(comp)| - the SCALE of a conserved
  // total, which is what a conservation error has to be measured against. The
  // signed total is the wrong yardstick for e: it is a sum over positive and
  // negative formation energies, so it can sit near zero by cancellation and
  // make a harmless absolute error look enormous
  amrex::Real sumMassWeightedAbs(int lev, int comp) const;

  // Fernando-Clem: per-cell signed face mass fluxes [g/s] consumed by the
  // splicing operator (component = advection::Face; positive = incoming,
  // negative = outgoing). Filled by the LES coupling each step
  void defineFluxContainer(int lev);
  amrex::MultiFab& fluxFaces() { return m_flux_faces; }
  const amrex::MultiFab& fluxFaces() const { return m_flux_faces; }

  // Fernando-Clem: physical BC types per direction (PCPhysBCType encoding:
  // inflow == 1, outflow == 2); needed by splicing at domain boundaries
  void setPhysBC(const amrex::BCRec& bc)
  {
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      m_lo_bc[dir] = bc.lo(dir);
      m_hi_bc[dir] = bc.hi(dir);
    }
  }
  const amrex::IntVect& loBC() const { return m_lo_bc; }
  const amrex::IntVect& hiBC() const { return m_hi_bc; }

private:
  std::unique_ptr<ClemParticleContainer> m_pc;
  // Fernando-Clem: one map per AMR level, keyed by grid index
  amrex::Vector<CellVectors> m_cell_vectors;
  bool m_cell_map_valid = false;

  // Fernando-Clem: face mass-flux buffer (2*AMREX_SPACEDIM comps per cell)
  amrex::MultiFab m_flux_faces;
  // Fernando-Clem: physical BC types (PCPhysBCType); 0 = interior default
  amrex::IntVect m_lo_bc{AMREX_D_DECL(0, 0, 0)};
  amrex::IntVect m_hi_bc{AMREX_D_DECL(0, 0, 0)};

  static int m_verbose;
  static bool m_freeze_species;
  static int m_plot_int;
  static bool m_pressure_coupling;
  static bool m_couple_back;
  static int m_couple_back_energy;
  static bool m_vol_renorm;
  static int m_couple_back_mode;
  static int m_regrid_int;
  static bool m_do_diffusion;
  static amrex::Real m_diffusion_cfl;
  static int m_diffusion_max_substeps;
  static bool m_do_react;
};

} // namespace clem
#endif
