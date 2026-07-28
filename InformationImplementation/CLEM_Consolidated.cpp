// ================================================================================
// CLEM (Compressible Linear-Eddy-Model LES) -- CONSOLIDATED SOURCE BUNDLE
// ================================================================================
// This is a READ-ONLY review artifact, NOT part of the build. It concatenates
// every file in Source/Clem/ (headers then implementations) so a fresh AI
// session/agent can load the FULL current implementation in one file instead
// of traversing ~10 files. Original paths are marked at each file boundary --
// edit the ORIGINALS there, never this bundle. Regenerate by re-concatenating.
//
// -------------------------------------------------------------------------
// WHAT THIS IS
// -------------------------------------------------------------------------
// CLEM-LES per Maxwell (2016 PhD thesis, "Turbulent Combustion Modelling of
// Fast-flames and Detonations Using Compressible LEM-LES"): each LES
// (supergrid) cell owns a 1-D Linear-Eddy-Model line of n_lem Lagrangian
// mass-coordinate elements ("particles" in AMReX terms) representing the
// sub-cell structure. Two-way coupling: LES->SGS pressure work (isentropic),
// SGS->LES filtered composition/energy write-back. Full docs:
// InformationImplementation/MaxwellCLEM_FaithfulBaseline.md (Maxwell's model,
// section-numbered [M x.xx]) and CouplingLESwithCLEM.md (multi-species
// generalization notes).
//
// SCOPE DECISION (binding, see MaxwellCLEM_FaithfulBaseline.md Sec 0/10):
// phase-1 is Maxwell-faithful with ONE accepted deviation -- multi-species
// real-gas chemistry (PelePhysics EOS + LiDryer mechanism) instead of
// Maxwell's calorically-perfect single-reactant gas. THIS DEVIATION IS THE
// ROOT OF EVERY OUTSTANDING PROBLEM BELOW: Maxwell's gas has T = p/rho
// identically, so matching pressure automatically matches internal energy --
// there is no separate channel for them to disagree. Real-gas rho*e =
// rho*e(p,T,Y) breaks that identity, so every composition write-back opens a
// gap between Y and rho*e that Maxwell's formulation never has to close.
//
// -------------------------------------------------------------------------
// CURRENT STATE (2026-07-27) -- READ THIS BEFORE CHANGING ANYTHING
// -------------------------------------------------------------------------
// Best-known-good CLEM config (baked into the canonical input files, see
// below): clem.freeze_species=1 (Maxwell: LES does NOT advect rho*Y_k),
// clem.pressure_coupling=1, clem.couple_back=1, clem.couple_back_energy=1
// (mode 1, conserved-delta -- the ONLY validated-safe energy mode),
// clem.vol_renorm=1 (ensemble-volume renormalization). do_diffusion=0 for
// pure advection / =1 when pelec.diffuse_*=1 to match.
//
// CANONICAL, SELF-DOCUMENTED INPUT FILES (use these, not ad-hoc configs):
//   CLEM side  Exec/RegTests/CLEM1DFlame/inputs-CLEM-ADV.inp   -> RESULTS/pltADV_CLEM_*
//   CLEM side  Exec/RegTests/CLEM1DFlame/inputs-CLEM-DIFF.inp  -> RESULTS/pltDIF_CLEM_*
//   pristine   ../1DFlame/inputs-REF-ADV.inp   (unmodified PeleC) -> ANALYSIS/pltADV_*
//   pristine   ../1DFlame/inputs-REF-DIFF.inp                     -> ANALYSIS/pltDIF_*
// (pristine tree = /home/fernando/FPCE/PeleC, separate checkout with NO Clem
// source at all -- the ground truth stock-PeleC reference.)
//
// MEASURED RESULTS (1000 steps, 256 cells, phi=0.4 H2/air LiDryer flame):
//   Advection-only:  velocity error +1.68 cm/s @500, +1.38 @1000 (bounded,
//                     small bulge at the flame front x~0.345 cm)
//   Advection+diffusion: velocity error -26.79 cm/s @500, -26.50 @1000
//                     (bounded, NOT growing, deficit vs reference's sharp
//                     viscous-heating spike at the front)
//   Coupling health in both: species-mass projection residual ~1e-15 to
//   8.9e-16 (machine precision, sum_k rho*Y_k = rho_LES enforced exactly);
//   zero NaN/abort/blow-up over 1000 steps.
//
// -------------------------------------------------------------------------
// THE OUTSTANDING PROBLEM -- PRECISE DIAGNOSIS
// -------------------------------------------------------------------------
// The residual velocity error (both cases) is NOT an energy-bookkeeping bug.
// It is a TRANSPORT-ORDER MISMATCH: splicing (SplicingOperator, below) moves
// parcels 1st-order (whole/split parcel transfer at exact face mass flux, no
// slope reconstruction) while the LES advects the same field with a 2nd-order
// Godunov/MUSCL scheme (Hydro.cpp, outside this bundle). From an identical
// baseline (measured offset ~6.6e-16, i.e. exactly synced) one CLEM step
// diverges from the matching LES step by ~2.3e-8 (pure advection, per-step,
// at the front) to ~1e-4 relative in Y(H2) after many steps -- this gap is
// IDENTICAL in the advection-only and the diffusion case (diffusion does not
// add to it; the LesCoupling::conservedConsistency instrument proved the
// species channel is otherwise clean, gap ~6e-12 in pure advection).
//
// The deficit's MECHANISM in the diffusion case: the write-back's small (~
// 1e-4) composition shift at the mismatched-transport EOS moves T by a few
// K; with pelec.diffuse_temp/enth=1 the LES CONDUCTION operator then acts on
// that shifted T every subsequent step, amplifying a transport artifact into
// a sustained (but bounded) velocity deficit. This is WHY it's ~15x larger
// with diffusion on (-26.5) than pure advection (+1.4): conduction is a
// feedback amplifier on the same root cause, not a second, independent bug.
//
// EXHAUSTED AVENUE -- energy-channel patches (all measured, see
// LesCoupling::writeBackProjected / writeBackToLes below for the code):
//   mode 0 (energy untouched):        -30.08 @500, -29.82 @1000 (bounded)
//   mode 1 (conserved-delta, +renorm): -26.79 @500, -26.50 @1000 (BEST, ~11%
//                                       better than mode 0, still bounded)
//   mode 2 ("keep-T" reconciliation):  +52.7 @500, +78.3 @1000 (GROWS -- the
//                                       formation-energy of dY injected at
//                                       fixed T each step is a non-
//                                       conservative source -> pumps)
//   pressure-neutral (earliest, code removed): killed the bulge locally but
//                                       is non-conservative -- injects an
//                                       energy pulse whose acoustic wake
//                                       grows LINEARLY. Also rejected.
// Conclusion: no energy formula can be smaller than the transport error that
// creates the need for it. Mode 1 is the ceiling of that search.
//
// -------------------------------------------------------------------------
// OPEN DIRECTIONS (not yet started; pick one, each is a real scope of work)
// -------------------------------------------------------------------------
//  1. CONVERGENCE STUDY (cheap, no new code): rerun inputs-CLEM/REF-DIFF at
//     512/1024 cells. If the -26.5 deficit shrinks with resolution as
//     expected for a 1st-order transport error, it proves this is a
//     controlled, convergent discretization artifact (same class as
//     Maxwell's own documented splicing-diffusion limitation #2, see
//     MaxwellCLEM_FaithfulBaseline.md Sec 9), not a modeling bug.
//  2. HIGHER-ORDER SPLICE (real fix, real work): rewrite
//     SplicingOperator::advectCells (below) with a slope-limited /
//     piecewise-linear parcel reconstruction so its transport order matches
//     the LES's Godunov scheme. Directly targets the measured root cause.
//  3. SIMPLIFIED-EOS ISOLATION TEST (evidence, not a fix): build a small
//     synthetic harness with a calorically-perfect two-state gas (T=p/rho
//     exactly) run through the SAME coupling machinery. If the deficit
//     vanishes there, it proves the residual is 100% attributable to the
//     multi-species EOS complexity, confirming the diagnosis with hard data
//     rather than inference.
//  4. IMPLEMENT REACTION (the actual point of CLEM): accept the current
//     bounded, fully-characterized deficit as the phase-1 baseline and wire
//     the subgrid reaction stage (currently a stub -- see the "(future
//     stages) stirring -> reaction" comment in Algorithm::advance below).
//     CLEM is presently a near no-op without reaction (diffusion-only
//     feedback measured ~0 by design -- the LEM line is PERIODIC per cell,
//     so diffusion conserves the filtered mean exactly, see ClemDiffusion.H
//     header comment). Reaction is where CLEM stops being a no-op and where
//     the whole coupling machinery above finally earns its keep.
//
// -------------------------------------------------------------------------
// FILE MANIFEST (concatenated below, in this reading order)
// -------------------------------------------------------------------------
//  ClemIndex.H            per-particle layout (RealData/IntData enums),
//                          compile-time n_lem and splice-policy config
//  ClemParticleContainer.H AMReX particle container typedefs (AoS storage)
//  ClemScratch.H           per-rank reusable line-operator scratch arena
//  ClemEosUtil.H           bias-free (rho,e,Y)->T EOS inversion helper
//  ClemAdvection.H         face/flux bookkeeping + periodic/BC helpers used
//                          by splicing
//  ClemManager.H/.cpp      owns the particle container + cell->particle map;
//                          all clem.* runtime params; conservation-sum
//                          helpers; particle plotfile I/O
//  ClemRegrid.H/.cpp       rebins each line back to n_lem EQUAL-MASS elements
//                          (Maxwell thesis p.172); T/Y mass-averaged, (e,p)
//                          EOS-RECONSTRUCTED (not averaged -- this was a
//                          fixed particle-blowup bug, see git history)
//  ClemSplicing.H/.cpp     INTER-CELL transport at the LES face mass fluxes
//                          (Maxwell steps 23-24). *** THE SUSPECT FILE for
//                          open direction #2 above -- advectCells() is the
//                          1st-order parcel-transfer scheme to replace ***
//  ClemDiffusion.H/.cpp    LINE-INTERNAL molecular diffusion (periodic,
//                          mass-coordinate, isobaric closure). Conserves the
//                          filtered mean exactly by construction -- inert to
//                          the LES until reaction exists (see header note)
//  ClemLesCoupling.H/.cpp  THE COUPLING OPERATORS -- isentropicPressureUpdate
//                          (LES->SGS), renormalizeEnsembleVolume (the vol_renorm
//                          fix), snapshotFilteredMean, writeBackProjected
//                          (Maxwell/freeze=1 path, THE energy-mode switch
//                          lives here), writeBackToLes (React.cpp-style
//                          increment path, freeze=0), conservedConsistency
//                          (the LES-vs-SGS diagnostic instrument that
//                          localized the transport-order gap)
//  ClemAlgorithm.H/.cpp    TOP-LEVEL ORCHESTRATION -- Algorithm::advance()
//                          composes every stage above in Maxwell's order;
//                          READ THIS FIRST if you want the one-function
//                          overview, then dive into whichever operator
// ================================================================================


// ================================================================================
// FILE: Source/Clem/ClemIndex.H  (original -- edit THERE, not in this bundle)
// ================================================================================
#ifndef CLEMINDEX_H
#define CLEMINDEX_H

#include <AMReX_REAL.H>
#include "mechanism.H"

// Fernando-Clem: compile-time number of LEM elements per LES cell
// (alpaca-style: override at build time with -DCLEM_NLEM=<n>)
#ifndef CLEM_NLEM
#define CLEM_NLEM 12
#endif

// Fernando-Clem: compile-time splicing policy for the strategy comparison
//   0 = side_aware (Maxwell steps 23-24: donor from the exit-facing end,
//       arrivals spliced onto the end matching their entry face)
//   1 = append_downstream (legacy: donor from the back, all arrivals appended
//       at the downstream end - known artificial front diffusivity)
#ifndef CLEM_SPLICE_POLICY
#define CLEM_SPLICE_POLICY 0
#endif

namespace clem {

namespace config {
// Fernando-Clem: compile-time configuration, accessed as constexpr everywhere
static constexpr int n_lem = CLEM_NLEM;

enum class SplicePolicy : int { side_aware = 0, append_downstream = 1 };
static constexpr SplicePolicy splice_policy =
  static_cast<SplicePolicy>(CLEM_SPLICE_POLICY);
} // namespace config

// Fernando-Clem: per-particle real components p_k = [rho, m, v, T, eint, p,
// Y_k] (primitive storage, Maxwell's element state). press is the element's
// own pressure: reference of the isentropic LES->SGS stage. xflux is
// transient splicing bookkeeping: the face mass flux that carried the
// parcel, used to order same-end arrivals; cleared on reindex
struct RealData
{
  static constexpr int rho = 0;   // density
  static constexpr int mass = 1;  // element mass
  static constexpr int vol = 2;   // element volume
  static constexpr int T = 3;     // temperature
  static constexpr int eint = 4;  // specific internal energy
  static constexpr int press = 5; // element pressure
  static constexpr int xflux = 6; // |face mass flux| of the last splice move
  static constexpr int Y0 = 7;    // first species mass fraction, Y_k = Y0 + k
  static constexpr int ncomps = 7 + NUM_SPECIES;
};

// Fernando-Clem: per-particle int components. entry_side is transient
// splicing bookkeeping (0 = resident, 1 = entered through the left end,
// 2 = entered right, 3 = transverse arrival); cleared on reindex
struct IntData
{
  static constexpr int lem_index = 0;  // position along the 1-D LEM line
  static constexpr int entry_side = 1; // splicing entry code, see above
  static constexpr int ncomps = 2;
};

} // namespace clem
#endif

// ================================================================================
// FILE: Source/Clem/ClemParticleContainer.H  (original -- edit THERE, not in this bundle)
// ================================================================================
#ifndef CLEMPARTICLECONTAINER_H
#define CLEMPARTICLECONTAINER_H

#include <AMReX_AmrParticles.H>

#include "ClemIndex.H"

namespace clem {

// Fernando-Clem: iterators mirror the spray pattern (MyParIter/MyParConstIter)
class ClemParIter
  : public amrex::ParIter<RealData::ncomps, IntData::ncomps, 0, 0>
{
public:
  using amrex::ParIter<RealData::ncomps, IntData::ncomps, 0, 0>::ParIter;
};

class ClemParConstIter
  : public amrex::ParConstIter<RealData::ncomps, IntData::ncomps, 0, 0>
{
public:
  using amrex::ParConstIter<RealData::ncomps, IntData::ncomps, 0, 0>::
    ParConstIter;
};

// Fernando-Clem: storage-only container (AoS, no SoA comps, like spray);
// physics operators live outside so strategies can be swapped at compile time
class ClemParticleContainer
  : public amrex::AmrParticleContainer<RealData::ncomps, IntData::ncomps, 0, 0>
{
public:
  using BaseType =
    amrex::AmrParticleContainer<RealData::ncomps, IntData::ncomps, 0, 0>;
  using ParticleType = BaseType::ParticleType;

  explicit ClemParticleContainer(amrex::AmrCore* amr) : BaseType(amr) {}

  // Fernando-Clem: single-level ctor so the container is testable without a
  // full Amr hierarchy
  ClemParticleContainer(
    const amrex::Geometry& geom,
    const amrex::DistributionMapping& dmap,
    const amrex::BoxArray& ba)
    : BaseType({geom}, {dmap}, {ba}, {})
  {
  }
};

} // namespace clem
#endif

// ================================================================================
// FILE: Source/Clem/ClemScratch.H  (original -- edit THERE, not in this bundle)
// ================================================================================
#ifndef CLEMSCRATCH_H
#define CLEMSCRATCH_H

#include <AMReX_Vector.H>

#include "ClemIndex.H"

namespace clem {

// Fernando-Clem: per-rank reusable scratch arena for the LEM line operators.
//
// A line operator cannot work on the particle AoS in place: it needs the cell's
// elements gathered in LINE order (the map hands them out in AoS order), it
// needs derived state that is not stored on the particle (X_k, h_k, rho*D_k,
// rho*lambda, Cp), and it needs face and rate buffers that have no particle
// counterpart at all. That is a lot of copied information per cell, and sizing
// it per cell would put n_cell * n_substep heap allocations inside the inner
// loop of every step.
//
// So the buffers are owned by ONE arena per rank, grown on the first cell of
// the run and reused unchanged for every cell of every step afterwards (resize
// never shrinks). The arena is thread_local, so if a line operator is ever put
// inside an OpenMP region each thread gets its own line and the buffers stay
// race-free by construction; with MPI-only there is exactly one arena per rank.
//
// Layout: element quantities are indexed [l], per-species quantities are
// species-major per element, [l * NUM_SPECIES + k]; face quantities are indexed
// by the interior face f, which sits between elements f and f+1 (so face f is
// the "l+1/2" face of element l = f). The LEM line is PERIODIC (Maxwell's
// homogeneous sub-grid sample): face n_elem-1 wraps element n_elem-1 -> 0, so
// there are n_elem faces. Inter-cell transport stays the splicing operator's job.
struct LineScratch
{
  // Fernando-Clem: the per-rank instance; all line operators share it
  static LineScratch& get()
  {
    static thread_local LineScratch arena;
    return arena;
  }

  // Fernando-Clem: gathered element state (line order, lem_index ascending)
  amrex::Vector<int> pidx;          // index of the element in the tile AoS
  amrex::Vector<amrex::Real> mass;  // element mass = the mass-coordinate dm
  amrex::Vector<amrex::Real> rho;
  amrex::Vector<amrex::Real> T;
  amrex::Vector<amrex::Real> eint;
  amrex::Vector<amrex::Real> press; // element pressure, held constant over the
                                    // sub-cycle (the LEM array is isobaric)
  amrex::Vector<amrex::Real> Y;

  // Fernando-Clem: derived state, recomputed from (rho, T, Y) every sub-step so
  // the driving forces and the coefficients stay consistent as the line evolves
  amrex::Vector<amrex::Real> X;      // mole fractions
  amrex::Vector<amrex::Real> hk;     // species enthalpies (formation included)
  amrex::Vector<amrex::Real> rhoD;   // rho * (rho D_k W_k/W) - the mass-coordinate
                                     // species coefficient
  amrex::Vector<amrex::Real> rholam; // rho * lambda - the mass-coordinate
                                     // thermal coefficient
  amrex::Vector<amrex::Real> cp;     // mixture Cp
  amrex::Vector<amrex::Real> cpk;    // per-species Cp_k; the temperature
                                     // equation's interspecies term uses Cp_k,
                                     // NOT h_k (that would be the energy form)
  amrex::Vector<amrex::Real> wbar;   // mean molecular weight; the species
                                     // stability limit needs W/W_k (see
                                     // DiffusionOperator::stableSubStep)

  // Fernando-Clem: interior-face fluxes
  amrex::Vector<amrex::Real> Sk;   // corrected species mass flux (mass equation)
  amrex::Vector<amrex::Real> qf;   // conductive heat flux (rho lambda) dT/dm
  amrex::Vector<amrex::Real> qint; // per-face interspecies product
                                   // sum_k Cp_k S_k dT/dm (temperature equation)

  // Fernando-Clem: rates of the current sub-step
  amrex::Vector<amrex::Real> dYdt;
  amrex::Vector<amrex::Real> dTdt; // T is the transported thermal variable

  int n_elem = 0; // elements on the gathered line
  int n_face = 0; // faces = n_elem (periodic line: face n_elem-1 wraps to 0)

  // Fernando-Clem: grow to hold a line of n elements. Steady state is reached
  // on the first cell, so every later call is just the two assignments
  void resize(const int n)
  {
    n_elem = n;
    n_face = n; // periodic: one face per element (the last wraps n-1 -> 0)

    if (static_cast<int>(mass.size()) >= n) {
      return; // Fernando-Clem: already big enough - the common case
    }

    const int ns = n * NUM_SPECIES;
    const int nf = n; // periodic line

    pidx.resize(n);
    mass.resize(n);
    rho.resize(n);
    T.resize(n);
    eint.resize(n);
    press.resize(n);
    cp.resize(n);
    wbar.resize(n);
    rholam.resize(n);
    dTdt.resize(n);

    Y.resize(ns);
    X.resize(ns);
    hk.resize(ns);
    cpk.resize(ns);
    rhoD.resize(ns);
    dYdt.resize(ns);

    qf.resize(nf);
    qint.resize(nf);
    Sk.resize(nf * NUM_SPECIES);
  }
};

} // namespace clem
#endif

// ================================================================================
// FILE: Source/Clem/ClemEosUtil.H  (original -- edit THERE, not in this bundle)
// ================================================================================
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

// ================================================================================
// FILE: Source/Clem/ClemAdvection.H  (original -- edit THERE, not in this bundle)
// ================================================================================
#ifndef CLEMADVECTION_H
#define CLEMADVECTION_H

#include <utility>

#include <AMReX.H>
#include <AMReX_IntVect.H>
#include <AMReX_Box.H>

// Fernando-Clem: face/flux bookkeeping and boundary helpers for the splicing
// operator, ported from the old clem_advection.H

namespace clem {

namespace advection {

enum class Face { Left = 0, Right = 1, Bottom = 2, Top = 3, South = 4, North = 5 };

static constexpr int num_faces = 2 * AMREX_SPACEDIM;

struct FluxDir
{
  amrex::Real flux;
  Face face;
};

constexpr amrex::IntVect face_offsets[] = {
  amrex::IntVect{AMREX_D_DECL(-1, 0, 0)}, // Left
  amrex::IntVect{AMREX_D_DECL(1, 0, 0)},  // Right
  amrex::IntVect{AMREX_D_DECL(0, -1, 0)}, // Bottom
  amrex::IntVect{AMREX_D_DECL(0, 1, 0)},  // Top
#if AMREX_SPACEDIM == 3
  amrex::IntVect{AMREX_D_DECL(0, 0, -1)}, // South
  amrex::IntVect{AMREX_D_DECL(0, 0, 1)},  // North
#endif
};

AMREX_FORCE_INLINE
amrex::IntVect
neighbor_iv(const amrex::IntVect& iv, Face face)
{
  return iv + face_offsets[static_cast<int>(face)];
}

} // namespace advection

namespace boundary {

// Fernando-Clem: wrap a destination cell through periodic boundaries;
// second member reports whether wrapping occurred
AMREX_FORCE_INLINE
std::pair<amrex::IntVect, bool>
wrap_iv(
  const amrex::IntVect& iv,
  const amrex::Box& domain,
  const amrex::IntVect& is_periodic)
{
  amrex::IntVect r = iv;
  if (domain.contains(r)) {
    return {r, false};
  }
  bool was_wrapped = false;
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    if (is_periodic[dir] != 0) {
      const int length = domain.length(dir);
      if (r[dir] < domain.smallEnd(dir)) {
        r[dir] += length;
        was_wrapped = true;
      } else if (r[dir] > domain.bigEnd(dir)) {
        r[dir] -= length;
        was_wrapped = true;
      }
    }
  }
  return {r, was_wrapped};
}

// Fernando-Clem: cell outside the domain across a non-periodic outflow
// boundary (PCPhysBCType encoding: outflow == 2)
AMREX_FORCE_INLINE
bool
is_outside_domain_on_outflow(
  const amrex::IntVect& iv,
  const amrex::Box& domain,
  const amrex::IntVect& is_periodic,
  const amrex::IntVect& lo_bc,
  const amrex::IntVect& hi_bc)
{
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    if (is_periodic[dir] != 0) {
      continue;
    }
    if (iv[dir] < domain.smallEnd(dir) && lo_bc[dir] == 2) {
      return true;
    }
    if (iv[dir] > domain.bigEnd(dir) && hi_bc[dir] == 2) {
      return true;
    }
  }
  return false;
}

// Fernando-Clem: crossing this face enters the domain through a non-periodic
// inflow boundary (PCPhysBCType encoding: inflow == 1); dst is the would-be
// destination (outside the domain for a boundary face)
AMREX_FORCE_INLINE
bool
is_inflow_boundary_face(
  const amrex::IntVect& dst,
  const amrex::Box& domain,
  const amrex::IntVect& is_periodic,
  const amrex::IntVect& lo_bc,
  const amrex::IntVect& hi_bc)
{
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    if (is_periodic[dir] != 0) {
      continue;
    }
    if (dst[dir] < domain.smallEnd(dir) && lo_bc[dir] == 1) {
      return true;
    }
    if (dst[dir] > domain.bigEnd(dir) && hi_bc[dir] == 1) {
      return true;
    }
  }
  return false;
}

} // namespace boundary
} // namespace clem
#endif

// ================================================================================
// FILE: Source/Clem/ClemManager.H  (original -- edit THERE, not in this bundle)
// ================================================================================
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
};

} // namespace clem
#endif

// ================================================================================
// FILE: Source/Clem/ClemManager.cpp  (original -- edit THERE, not in this bundle)
// ================================================================================
#include <cmath>

#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>

#include "PelePhysics.H"

#include "ClemManager.H"

namespace clem {

int ClemManager::m_verbose = 0;
bool ClemManager::m_freeze_species = false;
int ClemManager::m_plot_int = -1;
// Fernando-Clem: coupling stages default OFF so a run without clem.* inputs
// is transport-only and the LES stays bit-identical to stock PeleC
bool ClemManager::m_pressure_coupling = false;
bool ClemManager::m_couple_back = false;
int ClemManager::m_couple_back_energy = 0;
bool ClemManager::m_vol_renorm = true;
int ClemManager::m_regrid_int = 1;
bool ClemManager::m_do_diffusion = false;
// Fernando-Clem: 1/2 is the stability limit of the 2nd-order central stencil
amrex::Real ClemManager::m_diffusion_cfl = 0.5;
int ClemManager::m_diffusion_max_substeps = 1000;

void
ClemManager::readParams()
{
  amrex::ParmParse pp("clem");
  pp.query("verbose", m_verbose);
  pp.query("freeze_species", m_freeze_species);
  pp.query("plot_int", m_plot_int);
  pp.query("pressure_coupling", m_pressure_coupling);
  pp.query("couple_back", m_couple_back);
  pp.query("couple_back_energy", m_couple_back_energy);
  pp.query("vol_renorm", m_vol_renorm);
  // Fernando-Clem: both alternative substitution modes were measured and
  // rejected - 1 (isobaric) pumps formation enthalpy and blows up; 2
  // (constant-T) double-counts the chemical enthalpy the LES energy flux
  // already transports (T error 20x mode 0, linear growth). Fixed-e
  // substitution (0) is the conservative, correct coupling
  pp.query("regrid_int", m_regrid_int);

  pp.query("do_diffusion", m_do_diffusion);
  pp.query("diffusion_cfl", m_diffusion_cfl);
  pp.query("diffusion_max_substeps", m_diffusion_max_substeps);

  if (m_diffusion_cfl <= 0.0 || m_diffusion_cfl > 0.5) {
    amrex::Abort(
      "clem.diffusion_cfl must be in (0, 0.5]: 0.5 is the stability limit of "
      "the 2nd-order central stencil used on the LEM line.");
  }
  if (m_diffusion_max_substeps < 1) {
    amrex::Abort("clem.diffusion_max_substeps must be >= 1.");
  }

  if (m_verbose > 0) {
    amrex::Print() << "CLEM: n_lem = " << config::n_lem
                   << ", real comps = " << RealData::ncomps
                   << ", int comps = " << IntData::ncomps
                   << ", freeze_species = " << m_freeze_species
                   << ", pressure_coupling = " << m_pressure_coupling
                   << ", couple_back = " << m_couple_back
                   << ", do_diffusion = " << m_do_diffusion << '\n';
  }
}

void
ClemManager::define(amrex::AmrCore* amr)
{
  AMREX_ALWAYS_ASSERT(m_pc == nullptr);
  m_pc = std::make_unique<ClemParticleContainer>(amr);
  m_pc->SetVerbose(m_verbose);
}

void
ClemManager::defineSingleLevel(
  const amrex::Geometry& geom,
  const amrex::DistributionMapping& dmap,
  const amrex::BoxArray& ba)
{
  AMREX_ALWAYS_ASSERT(m_pc == nullptr);
  m_pc = std::make_unique<ClemParticleContainer>(geom, dmap, ba);
  m_pc->SetVerbose(m_verbose);
}

void
ClemManager::initParticles(const int lev)
{
  BL_PROFILE("ClemManager::initParticles()");
  AMREX_ALWAYS_ASSERT(m_pc != nullptr);

  const auto& geom = m_pc->Geom(lev);
  const auto* dx = geom.CellSize();
  const auto* plo = geom.ProbLo();
  // Fernando-Clem: LEM elements evenly spaced along x inside each cell so the
  // cell->index relation is unambiguous
  const amrex::Real dx_inner =  dx[0] / static_cast<amrex::Real>(config::n_lem + 1);

  for (amrex::MFIter mfi = m_pc->MakeMFIter(lev, false); mfi.isValid(); ++mfi) {
    const amrex::Box& bx  = mfi.validbox();
    const int grid_id     = mfi.index();
    const int tile_id     = mfi.LocalTileIndex();
    auto& ptile           = m_pc->GetParticles(lev)[{grid_id, tile_id}];

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      for (int ip = 0; ip < config::n_lem; ++ip) {
        ClemParticleContainer::ParticleType p;
        p.id()  = ClemParticleContainer::ParticleType::NextID();
        p.cpu() = amrex::ParallelDescriptor::MyProc();

        AMREX_D_TERM(
            p.pos(0) =   plo[0] + iv[0] * dx[0] + static_cast<amrex::Real>(ip + 1) * dx_inner;
          , p.pos(1) =  plo[1] + (iv[1] + 0.5) * dx[1];
          , p.pos(2) =  plo[2] + (iv[2] + 0.5) * dx[2];)
        // Fernando-Clem: placeholder state; the coupling fills
        // [rho, m, v, T, eint, Y_k] from the LES cell
        for (int n = 0; n < RealData::ncomps; ++n) {
          p.rdata(n) = 0.0;
        }
        p.idata(IntData::lem_index) = ip;
        p.idata(IntData::entry_side) = 0;
        ptile.push_back(p);
      }
    }
  }
  m_cell_map_valid = false;
}

void
ClemManager::setParticlesFromState(
  const int lev,
  const amrex::MultiFab& state,
  const int rho_indx,
  const int temp_indx,
  const int eint_indx,
  const int spec_indx)
{
  BL_PROFILE("ClemManager::setParticlesFromState()");
  AMREX_ALWAYS_ASSERT(m_pc != nullptr);
  // Fernando-Clem: state and particles must live on the same grids
  AMREX_ALWAYS_ASSERT(state.boxArray() == m_pc->ParticleBoxArray(lev));

  const auto* dx = m_pc->Geom(lev).CellSize();
  amrex::Real cell_vol = 1.0;
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    cell_vol *= dx[d];
  }
  // Fernando-Clem: equal volume share per element so sum_p m_p = rho * V_cell
  const amrex::Real vol_p = cell_vol / static_cast<amrex::Real>(config::n_lem);

  auto eos = pele::physics::PhysicsType::eos();

  for (ClemParIter pti(*m_pc, lev); pti.isValid(); ++pti) {
    const auto& sarr = state.const_array(pti.index());
    auto& particles  = pti.GetArrayOfStructs();
    const int np     = pti.numParticles();

    for (int pindex = 0; pindex < np; ++pindex) {
      auto& p                 = particles[pindex];
      const amrex::IntVect iv = m_pc->Index(p, lev);
      const amrex::Real rho   = sarr(iv, rho_indx);

      p.rdata(RealData::rho)  = rho;
      p.rdata(RealData::vol)  = vol_p;
      p.rdata(RealData::mass) = rho * vol_p;
      p.rdata(RealData::T)    = sarr(iv, temp_indx);
      // Fernando-Clem: state stores rho*e and rho*Y_k -> convert to specific
      p.rdata(RealData::eint) = sarr(iv, eint_indx) / rho;
      amrex::Real Y[NUM_SPECIES];
      for (int k = 0; k < NUM_SPECIES; ++k) {
        Y[k] = sarr(iv, spec_indx + k) / rho;
        p.rdata(RealData::Y0 + k) = Y[k];
      }
      // Fernando-Clem: element pressure = EOS pressure of its own state
      // (equals the cell pressure at t=0)
      amrex::Real P = 0.0;
      eos.RTY2P(rho, p.rdata(RealData::T), Y, P);
      p.rdata(RealData::press) = P;
    }
  }
}

void
ClemManager::rebuildCellMap(const int lev)
{
  BL_PROFILE("ClemManager::rebuildCellMap()");
  AMREX_ALWAYS_ASSERT(m_pc != nullptr);

  if (m_cell_vectors.size() < lev + 1) {
    m_cell_vectors.resize(lev + 1);
  }

  // Fernando-Clem: reset ownership lists on the valid boxes of this level
  for (amrex::MFIter mfi = m_pc->MakeMFIter(lev, false); mfi.isValid(); ++mfi) {
    auto& fab = m_cell_vectors[lev][mfi.index()];
    fab.resize(mfi.validbox());
    fab.setVal(std::vector<int>{});
  }

  // Fernando-Clem: bin each particle into its cell (index into the grid AoS).
  // Untiled layout assumed so the AoS index is grid-global
  for (ClemParIter pti(*m_pc, lev); pti.isValid(); ++pti) {
    AMREX_ALWAYS_ASSERT(pti.LocalTileIndex() == 0);
    const int np      = pti.numParticles();
    const int grid_id = pti.index();
    auto& particles   = pti.GetArrayOfStructs();
    auto& fab         = m_cell_vectors[lev][grid_id];

    for (int pindex = 0; pindex < np; ++pindex) {
      const auto& p           = particles[pindex];
      const amrex::IntVect iv = m_pc->Index(p, lev);
      fab(iv).push_back(pindex);
    }
  }
  m_cell_map_valid = true;
}

const std::vector<int>&
ClemManager::particlesInCell(
  const int lev, const int grid_id, const amrex::IntVect& iv) const
{
  AMREX_ALWAYS_ASSERT(m_cell_map_valid);
  return m_cell_vectors[lev].at(grid_id)(iv);
}

bool
ClemManager::verifyCellMap(const int lev) const
{
  BL_PROFILE("ClemManager::verifyCellMap()");
  AMREX_ALWAYS_ASSERT(m_pc != nullptr);

  bool ok = m_cell_map_valid;
  amrex::Long n_mapped = 0;

  if (ok) {
    for (amrex::MFIter mfi = m_pc->MakeMFIter(lev, false); mfi.isValid();++mfi) {
      const amrex::Box& bx  = mfi.validbox();
      const auto& fab       = m_cell_vectors[lev].at(mfi.index());

      for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
        const auto ncell = static_cast<int>(fab(iv).size());
        n_mapped += ncell;
        // Fernando-Clem: structural invariant - exactly n_lem elements per cell
        if (ncell != config::n_lem) {
          ok = false;
          amrex::AllPrint()
            << "CLEM: cell " << iv << " on grid " << mfi.index() << " owns "
            << ncell << " particles, expected " << config::n_lem << '\n';
        }
      }
    }
    // Fernando-Clem: no particle lost or double-counted on this rank
    const amrex::Long n_local = m_pc->NumberOfParticlesAtLevel(lev, true, true);
    if (n_mapped != n_local) {
      ok = false;
      amrex::AllPrint() << "CLEM: mapped " << n_mapped << " particles but "
                        << n_local << " are stored on this rank\n";
    }
  }

  amrex::ParallelDescriptor::ReduceBoolAnd(ok);
  if (m_verbose > 0) {
    amrex::Print() << "CLEM: cell map verification "
                   << (ok ? "passed" : "FAILED") << " on level " << lev
                   << '\n';
  }
  return ok;
}

void
ClemManager::redistribute(const int lbase)
{
  AMREX_ALWAYS_ASSERT(m_pc != nullptr);
  m_pc->Redistribute(lbase);
  m_cell_map_valid = false;
}

void
ClemManager::writePlotFile(
  const std::string& dir, const std::string& name) const
{
  AMREX_ALWAYS_ASSERT(m_pc != nullptr);
  const int lev = 0;

  // Fernando-Clem: ParaView cannot read a bare particle directory; it needs a
  // grid plotfile that contains the particle data as a subdirectory (same
  // layout spray uses). Write particle counts as the cell field, then embed
  // the particles under <dir>/<name>.
  amrex::MultiFab cnt(
    m_pc->ParticleBoxArray(lev), m_pc->ParticleDistributionMap(lev), 1, 0);
  cnt.setVal(0.0);
  m_pc->Increment(cnt, lev);
  amrex::WriteSingleLevelPlotfile(
    dir, cnt, {"particle_count"}, m_pc->Geom(lev), 0.0, 0);

  amrex::Vector<std::string> real_names = {
    "density", "mass", "volume", "temperature", "eint", "pressure",
    "mass_flux"};
  amrex::Vector<std::string> spec_names;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(spec_names);
  for (int n = 0; n < NUM_SPECIES; ++n) {
    real_names.push_back("Y_" + spec_names[n]);
  }
  const amrex::Vector<std::string> int_names = {"lem_index", "entry_side"};
  m_pc->WritePlotFile(dir, name, real_names, int_names);
}

amrex::Long
ClemManager::totalParticleCount(const int lev, const bool only_local) const
{
  AMREX_ALWAYS_ASSERT(m_pc != nullptr);
  return m_pc->NumberOfParticlesAtLevel(lev, true, only_local);
}

amrex::Real
ClemManager::sumParticleReal(const int lev, const int comp) const
{
  AMREX_ALWAYS_ASSERT(m_pc != nullptr);
  amrex::Real total = 0.0;
  for (ClemParConstIter pti(*m_pc, lev); pti.isValid(); ++pti) {
    const auto& particles = pti.GetArrayOfStructs();
    const int np          = pti.numParticles();
    for (int pindex = 0; pindex < np; ++pindex) {
      total += particles[pindex].rdata(comp);
    }
  }
  amrex::ParallelDescriptor::ReduceRealSum(total);
  return total;
}

void
ClemManager::defineFluxContainer(const int lev)
{
  AMREX_ALWAYS_ASSERT(m_pc != nullptr);
  // Fernando-Clem: one signed mass flux per face of every cell; the LES
  // coupling fills it, the splicing operator consumes it
  m_flux_faces.define(
    m_pc->ParticleBoxArray(lev), m_pc->ParticleDistributionMap(lev),
    2 * AMREX_SPACEDIM, 0);
  m_flux_faces.setVal(0.0);
}

amrex::Real
ClemManager::maxCellVolumeDefect(const int lev) const
{
  AMREX_ALWAYS_ASSERT(m_pc != nullptr);
  AMREX_ALWAYS_ASSERT(m_cell_map_valid);

  const auto* dx = m_pc->Geom(lev).CellSize();
  amrex::Real cell_vol = 1.0;
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    cell_vol *= dx[d];
  }
  const amrex::Real one_vol = 1.0 / cell_vol;

  amrex::Real defect = 0.0;
  for (ClemParConstIter pti(*m_pc, lev); pti.isValid(); ++pti) {
    const int grid_id = pti.index();
    const amrex::Box& bx = pti.validbox();
    const auto& particles = pti.GetArrayOfStructs();
    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      const auto& plist = particlesInCell(lev, grid_id, iv);
      amrex::Real v_sum = 0.0;
      for (const int pidx : plist) {
        v_sum += particles[pidx].rdata(RealData::vol);
      }
      defect = amrex::max<amrex::Real>(
        defect, std::abs(v_sum * one_vol - 1.0));
    }
  }
  amrex::ParallelDescriptor::ReduceRealMax(defect);
  return defect;
}

amrex::Real
ClemManager::sumMassWeighted(const int lev, const int comp) const
{
  AMREX_ALWAYS_ASSERT(m_pc != nullptr);
  amrex::Real total = 0.0;

  for (ClemParConstIter pti(*m_pc, lev); pti.isValid(); ++pti) {
    const auto& particles = pti.GetArrayOfStructs();
    const int np          = pti.numParticles();
    for (int pindex = 0; pindex < np; ++pindex) {
      total += particles[pindex].rdata(RealData::mass) *  particles[pindex].rdata(comp);
    }
  }
  amrex::ParallelDescriptor::ReduceRealSum(total);
  return total;
}

amrex::Real
ClemManager::sumMassWeightedAbs(const int lev, const int comp) const
{
  AMREX_ALWAYS_ASSERT(m_pc != nullptr);
  amrex::Real total = 0.0;
  for (ClemParConstIter pti(*m_pc, lev); pti.isValid(); ++pti) {
    const auto& particles   = pti.GetArrayOfStructs();
    const int np            = pti.numParticles();
    for (int pindex = 0; pindex < np; ++pindex) {
      total += particles[pindex].rdata(RealData::mass) * std::abs(particles[pindex].rdata(comp));
    }
  }
  amrex::ParallelDescriptor::ReduceRealSum(total);
  return total;
}

} // namespace clem

// ================================================================================
// FILE: Source/Clem/ClemRegrid.H  (original -- edit THERE, not in this bundle)
// ================================================================================
#ifndef CLEMREGRID_H
#define CLEMREGRID_H

#include "ClemManager.H"

namespace clem {

// Fernando-Clem: regridding operator (Maxwell thesis, p. 172). Rebins each
// LES cell's LEM line into exactly config::n_lem equal-mass elements while
// conserving mass, volume, species mass and internal energy. Stateless
// operator class (works through the manager's public API) so alternative
// regridding strategies can be swapped at compile time.
class Regridder
{
public:
  static void regrid(ClemManager& mgr, int lev);

private:
  // Fernando-Clem: top up cells owning fewer than n_lem elements with
  // zero-mass placeholders; they become merge targets and inherit mass
  static void ensureCellCapacity(ClemManager& mgr, int lev);

  // Fernando-Clem: per-cell two-pointer source->target equal-mass merge.
  // Note the particle state is primitive (T, specific e, mass fractions Y_k),
  // so the merge is a direct mass-weighted average - none of the old
  // conservative (rho Y, rho e) conversions are needed
  static void mergeCells(ClemManager& mgr, int lev);
};

} // namespace clem
#endif

// ================================================================================
// FILE: Source/Clem/ClemRegrid.cpp  (original -- edit THERE, not in this bundle)
// ================================================================================
#include <algorithm>

#include <AMReX_ParallelDescriptor.H>

#include "ClemEosUtil.H"
#include "ClemRegrid.H"

namespace clem {

void
Regridder::regrid(ClemManager& mgr, const int lev)
{
  BL_PROFILE("clem::Regridder::regrid()");
  AMREX_ALWAYS_ASSERT(mgr.isDefined());

  if (!mgr.IsCellMapValid()) {
    mgr.rebuildCellMap(lev);
  }
  ensureCellCapacity(mgr, lev);
  // Fernando-Clem: placeholders were pushed directly into the owning tile, so
  // no Redistribute is needed here - only the map must see them
  mgr.rebuildCellMap(lev);

  mergeCells(mgr, lev);

  // Fernando-Clem: Redistribute removes the id < 0 leftovers of the merge
  mgr.redistribute(lev);
  mgr.rebuildCellMap(lev);
}

void
Regridder::ensureCellCapacity(ClemManager& mgr, const int lev)
{
  BL_PROFILE("clem::Regridder::ensureCellCapacity()");
  auto& pc = mgr.particleContainer();
  const auto& geom = pc.Geom(lev);
  const auto* dx = geom.CellSize();
  const auto* plo = geom.ProbLo();
  const amrex::Real dx_inner = dx[0] / static_cast<amrex::Real>(config::n_lem + 1);

  for (amrex::MFIter mfi = pc.MakeMFIter(lev, false); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.validbox();
    const int grid_id = mfi.index();
    const int tile_id = mfi.LocalTileIndex();
    auto& ptile = pc.GetParticles(lev)[{grid_id, tile_id}];

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      // Fernando-Clem: only counts are read here (no AoS references are held
      // across push_back, which may reallocate the tile storage)
      const int n_in_cell = static_cast<int>(mgr.particlesInCell(lev, grid_id, iv).size());

      //if few particles, allocates aprticle tot have NUM_LEM
      for (int ip = n_in_cell; ip < config::n_lem; ++ip) {
        ClemParticleContainer::ParticleType p;
        p.id() = ClemParticleContainer::ParticleType::NextID();
        p.cpu() = amrex::ParallelDescriptor::MyProc();
        AMREX_D_TERM(
            p.pos(0)  = plo[0] + iv[0] * dx[0] + static_cast<amrex::Real>(ip + 1) * dx_inner;
          , p.pos(1)  = plo[1] + (iv[1] + 0.5) * dx[1];
          , p.pos(2)  = plo[2] + (iv[2] + 0.5) * dx[2];)

        for (int n = 0; n < RealData::ncomps; ++n) {
          p.rdata(n) = 0.0;
        }
        // Fernando-Clem: placeholders go to the end of the line so they land
        // in the tail target slots of the merge
        p.idata(IntData::lem_index) = ip;
        p.idata(IntData::entry_side) = 0;
        ptile.push_back(p);
      }
    }
  }
}

/**
 * @brief HEre the regridding process is implements, such as \Delta_m is enforced for 
 * the comming process 
 */
void
Regridder::mergeCells(ClemManager& mgr, const int lev)
{
  BL_PROFILE("clem::Regridder::mergeCells()");
  auto& pc                    = mgr.particleContainer();
  const auto& geom            = pc.Geom(lev);
  const auto* dx              = geom.CellSize();
  const auto* plo             = geom.ProbLo();
  const amrex::Real dx_inner  =   dx[0] / static_cast<amrex::Real>(config::n_lem + 1);

  constexpr int ncomps = RealData::ncomps;
  auto eos = pele::physics::PhysicsType::eos();
  // Fernando-Clem: scratch buffers reused across cells to avoid per-cell heap
  // allocations
  amrex::Vector<amrex::Real> src_data;
  std::vector<int> plist;

  for (ClemParIter pti(pc, lev); pti.isValid(); ++pti) {
    const int grid_id     = pti.index();
    const amrex::Box& bx  = pti.validbox();
    auto& particles       = pti.GetArrayOfStructs();

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      // Fernando-Clem: local sorted copy; the map's own lists stay untouched
      const auto& cell_list = mgr.particlesInCell(lev, grid_id, iv);
      plist.assign(cell_list.begin(), cell_list.end());
      std::sort(plist.begin(), plist.end(), [&](int a, int b) {
        return particles[a].idata(IntData::lem_index) <
               particles[b].idata(IntData::lem_index);
      });

      const int n_src = static_cast<int>(plist.size());
      AMREX_ALWAYS_ASSERT(n_src >= config::n_lem);

      // Fernando-Clem: snapshot the (primitive) source states; the same slots
      // are reused as merge targets below
      src_data.resize(n_src * ncomps);
      amrex::Real total_mass = 0.0;
      for (int i = 0; i < n_src; ++i) {
        const auto& p = particles[plist[i]];

        for (int n = 0; n < ncomps; ++n) {
          src_data[i * ncomps + n] = p.rdata(n);
        }
        total_mass += p.rdata(RealData::mass);
      }

      //delete container which total mass is lower than 0.0 or 0.0(very unlikely scenarion)
      if (total_mass <= 0.0) {
        // Fernando-Clem: empty cell (all placeholders) - nothing to merge;
        // keep n_lem placeholders, drop the rest
        for (int i = config::n_lem; i < n_src; ++i) {
          auto& p = particles[plist[i]];
          p.id() = -p.id();
        }
        continue;
      }

      const amrex::Real new_dm = total_mass / static_cast<amrex::Real>(config::n_lem);

      // Fernando-Clem: reset all slots; targets accumulate mass-weighted sums
      for (const int pidx : plist) {
        auto& p = particles[pidx];
        for (int n = 0; n < ncomps; ++n) {
          p.rdata(n) = 0.0;
        }
      }

      // Fernando-Clem: two-pointer walk - drain sources in line order into
      // equal-mass targets. T and Y_k accumulate as m*phi (mass-averaged in the
      // finalize); internal energy is NOT accumulated - it is reconstructed from
      // the merged (rho, T, Y). Volume transfers as mass/rho_src, conserved
      int is = 0;
      int it = 0;
      while (is < n_src && it < config::n_lem) {
        auto& p_t = particles[plist[it]];
        amrex::Real& m_src          = src_data[is * ncomps + RealData::mass];
        const amrex::Real rho_src   = src_data[is * ncomps + RealData::rho];
        const amrex::Real mass_need = new_dm - p_t.rdata(RealData::mass);

        const bool close_target     = (m_src >= mass_need && mass_need > 0.0);
        const amrex::Real dm_take   = close_target ? mass_need : m_src;
        // Fernando-Clem: volume of the chunk = mass / density of its source
        // element; accumulated so V_j = sum_i m_i/rho_i (Maxwell A.3). rho_src
        // is always > 0 here: real sources are drained before the loop can
        // reach the zero-mass placeholders (all targets fill first)
        const amrex::Real dv_take   = dm_take / rho_src;

        p_t.rdata(RealData::mass)   += dm_take;
        p_t.rdata(RealData::vol)    += dv_take;
        p_t.rdata(RealData::T)      += src_data[is * ncomps + RealData::T] * dm_take;
        // Fernando-Clem: press is NOT mass-averaged here - it is reconstructed
        // EOS-consistently from the merged (rho, T, Y) in the finalize below,
        // because it is the reference p_1 the next isentropic LES->SGS update
        // compresses from (a mass-weighted sum here -> tiny p_1 -> runaway).
        for (int sp = 0; sp < NUM_SPECIES; ++sp) {
          p_t.rdata(RealData::Y0 + sp) += src_data[is * ncomps + RealData::Y0 + sp] * dm_take;
        }

        m_src -= dm_take;
        if (close_target) {
          ++it;
        } else {
          ++is;
        }
      }

      // Fernando-Clem: all source mass must be consumed (FP dust tolerated)
      for (int i = 0; i < n_src; ++i) {
        AMREX_ASSERT_WITH_MESSAGE(
          std::abs(src_data[i * ncomps + RealData::mass]) < 1e-12 * new_dm,
          "CLEM regrid: source mass not fully consumed");
      }

      // Fernando-Clem: finalize targets - T and Y are mass-averaged (accumulated
      // m*phi divided by the new mass), rho from conserved mass/volume. T is the
      // transported/primary thermal variable, so internal energy is NOT averaged
      // but RECONSTRUCTED here from the merged (rho, T, Y): this keeps the regrid
      // T-consistent with the diffusion stage and lets the species formation
      // energies follow the merged composition. Trade-off: mass-averaging T does
      // not conserve sum(m e) exactly through the merge for a real gas
      for (int i = 0; i < config::n_lem; ++i) {
        auto& p = particles[plist[i]];
        const amrex::Real m = p.rdata(RealData::mass);
        if (m <= 0.0) {
          continue;
        }
        const amrex::Real one_m = 1.0 / m;
        
        p.rdata(RealData::T)                *= one_m;
        amrex::Real massfrac[NUM_SPECIES]    = {0.0};

        for (int sp = 0; sp < NUM_SPECIES; ++sp) {
          p.rdata(RealData::Y0 + sp) *= one_m;
          massfrac[sp] = p.rdata(RealData::Y0 + sp);
        }
        AMREX_ASSERT(p.rdata(RealData::vol) > 0.0);
        const amrex::Real rho = m / p.rdata(RealData::vol);
        p.rdata(RealData::rho) = rho;

        // Fernando-Clem: e AND press follow from the merged (rho, T, Y) by direct
        // EOS evaluation - NOT a mass-average. press MUST be EOS-consistent with
        // (rho, T, Y): the constant-pressure diffusion reconstructs rho from it
        // (PYT2RE), so a mass-averaged press injects a spurious density jump at
        // dt->0 and runs the ensemble volume away (the NaN we hit). It is also
        // the p_1 the next isentropic LES->SGS update compresses from
        amrex::Real e = 0.0;
        amrex::Real P = 0.0;
        eos.RTY2E(rho, p.rdata(RealData::T), massfrac, e);
        eos.RTY2P(rho, p.rdata(RealData::T), massfrac, P);
        p.rdata(RealData::eint)  = e;
        p.rdata(RealData::press) = P;
      }

      // Fernando-Clem: drop the drained slots beyond n_lem and re-index the
      // surviving line with evenly spaced positions
      for (int i = config::n_lem; i < n_src; ++i) {
        auto& p = particles[plist[i]];
        p.id() = -p.id();
      }
      for (int i = 0; i < config::n_lem; ++i) {
        auto& p = particles[plist[i]];
        AMREX_D_TERM(
            p.pos(0) = plo[0] + iv[0] * dx[0] + static_cast<amrex::Real>(i + 1) * dx_inner;
          , p.pos(1) = plo[1] + (iv[1] + 0.5) * dx[1];
          , p.pos(2) = plo[2] + (iv[2] + 0.5) * dx[2];)

        p.idata(IntData::lem_index) = i;
        AMREX_ASSERT_WITH_MESSAGE(std::abs(p.rdata(RealData::mass) - new_dm) < 1e-12 * new_dm, "CLEM regrid: element mass differs from the equal-mass target");
      }
    } // iv
  }   // pti
  mgr.invalidateCellMap();
}

} // namespace clem

// ================================================================================
// FILE: Source/Clem/ClemSplicing.H  (original -- edit THERE, not in this bundle)
// ================================================================================
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

// ================================================================================
// FILE: Source/Clem/ClemSplicing.cpp  (original -- edit THERE, not in this bundle)
// ================================================================================
#include <algorithm>
#include <vector>

#include <AMReX_ParallelDescriptor.H>

#include "PelePhysics.H"

#include "ClemAdvection.H"
#include "ClemSplicing.H"

namespace clem {

namespace {

using ParticleType = ClemParticleContainer::ParticleType;

// Fernando-Clem: donor end and receiver entry code for a given exit face,
// resolved at compile time from config::splice_policy
AMREX_FORCE_INLINE
bool
takeFromFront(const advection::Face face)
{
  if constexpr (config::splice_policy == config::SplicePolicy::side_aware) {
    return face == advection::Face::Left;
  } else {
    return false; // legacy: always drain the back of the line
  }
}

AMREX_FORCE_INLINE
int
entryCode(const advection::Face face)
{
  if constexpr (config::splice_policy == config::SplicePolicy::side_aware) {
    // material leaving through the donor's RIGHT face enters the receiver's
    // LEFT end (1); through the LEFT face -> RIGHT end (2); transverse -> 3
    return (face == advection::Face::Right) ? 1
           : (face == advection::Face::Left) ? 2
                                             : 3;
  } else {
    return 1; // legacy: every arrival appended the same way
  }
}

// Fernando-Clem: split-off parcel carries the (primitive) source state with
// the exported mass; volume from mass/rho keeps rho intensive
AMREX_FORCE_INLINE
void
makeSplitParticle(
  const amrex::IntVect& dst,
  const amrex::Real* plo,
  const amrex::Real* dx,
  const amrex::Real mass_out,
  const amrex::Real net_flux,
  const int entry_code,
  const ParticleType& p_in,
  ParticleType& p_out)
{
  AMREX_D_TERM(p_out.pos(0) = plo[0] + (dst[0] + 0.5) * dx[0];
               , p_out.pos(1) = plo[1] + (dst[1] + 0.5) * dx[1];
               , p_out.pos(2) = plo[2] + (dst[2] + 0.5) * dx[2];)

  for (int n = 0; n < RealData::ncomps; ++n) {
    p_out.rdata(n) = p_in.rdata(n);
  }
  p_out.rdata(RealData::mass) = mass_out;
  p_out.rdata(RealData::vol) = mass_out / p_in.rdata(RealData::rho);
  p_out.rdata(RealData::xflux) = std::abs(net_flux);

  // Fernando-Clem: reassigned by reindexAfterSplicing; large value keeps
  // arrivals behind residents within their entry group
  p_out.idata(IntData::lem_index) = 100000;
  p_out.idata(IntData::entry_side) = entry_code;
}

} // namespace

SpliceBalance
SplicingOperator::splice(
  ClemManager& mgr,
  const int lev,
  const amrex::Real dt,
  const amrex::MultiFab& state_bc,
  const int rho_indx,
  const int temp_indx,
  const int eint_indx,
  const int spec_indx)
{
  BL_PROFILE("clem::SplicingOperator::splice()");
  AMREX_ALWAYS_ASSERT(mgr.isDefined());
  AMREX_ALWAYS_ASSERT(mgr.fluxFaces().nComp() == advection::num_faces);
  AMREX_ALWAYS_ASSERT(state_bc.nGrow() >= 1);

  if (!mgr.IsCellMapValid()) {
    mgr.rebuildCellMap(lev);
  }

  SpliceBalance bal = advectCells(
    mgr, lev, dt, state_bc, rho_indx, temp_indx, eint_indx, spec_indx);

  // Fernando-Clem: moved parcels changed position (possibly rank); deleted
  // ones carry id < 0 and are purged here
  mgr.redistribute(lev);
  mgr.rebuildCellMap(lev);
  
  reindexAfterSplicing(mgr, lev);

  amrex::ParallelDescriptor::ReduceRealSum(bal.mass_in);
  amrex::ParallelDescriptor::ReduceRealSum(bal.mass_out);
  return bal;
}

SpliceBalance
SplicingOperator::advectCells(
  ClemManager& mgr,
  const int lev,
  const amrex::Real dt,
  const amrex::MultiFab& state_bc,
  const int rho_indx,
  const int temp_indx,
  const int eint_indx,
  const int spec_indx)
{
  auto& pc = mgr.particleContainer();
  const amrex::Geometry& geom     = pc.Geom(lev);
  const amrex::Real* plo          = geom.ProbLo();
  const amrex::Real* dx           = geom.CellSize();
  const amrex::Box& domain        = geom.Domain();

  AMREX_ALWAYS_ASSERT(state_bc.boxArray() == pc.ParticleBoxArray(lev));

  const amrex::IntVect is_periodic{AMREX_D_DECL(
    geom.isPeriodic(0) ? 1 : 0, geom.isPeriodic(1) ? 1 : 0,
    geom.isPeriodic(2) ? 1 : 0)};
  const amrex::IntVect& lo_bc = mgr.loBC();
  const amrex::IntVect& hi_bc = mgr.hiBC();

  SpliceBalance bal;
  std::vector<int> plist;

  for (ClemParIter pti(pc, lev); pti.isValid(); ++pti) {
    const int grid_id         = pti.index();
    const int tile_id         = pti.LocalTileIndex();
    const amrex::Box& bx      = pti.validbox();
    const auto& flux_arr      = mgr.fluxFaces().const_array(grid_id);
    const auto& state_arr     = state_bc.const_array(grid_id);
    auto& particle_tile       = pc.GetParticles(lev)[{grid_id, tile_id}];
    auto& particles           = particle_tile.GetArrayOfStructs();

    // Fernando-Clem: arrivals are appended after the cell loop; pushing here
    // would invalidate the AoS references
    amrex::Vector<ParticleType> new_particles;

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      // Fernando-Clem: face mass crossing over dt; sorted so the largest
      // export is served first when several faces drain the same line
      advection::FluxDir fluxes[advection::num_faces];
      for (int f = 0; f < advection::num_faces; ++f) {
        fluxes[f] = {flux_arr(iv, f) * dt, static_cast<advection::Face>(f)};
      }
      std::sort( std::begin(fluxes), std::end(fluxes), [](const advection::FluxDir& a, const advection::FluxDir& b) {
          return a.flux > b.flux;
        });

      // Fernando-Clem: local working copy in line order; the map's own lists
      // stay untouched (they are rebuilt after Redistribute)
      const auto& cell_list = mgr.particlesInCell(lev, grid_id, iv);
      plist.assign(cell_list.begin(), cell_list.end());

      std::sort(plist.begin(), plist.end(), [&](int a, int b) {
        return particles[a].idata(IntData::lem_index) <  particles[b].idata(IntData::lem_index);
      });

      for (int iface = 0; iface < advection::num_faces; ++iface) {
        const advection::Face face    = fluxes[iface].face;
        const amrex::Real net_flux    = fluxes[iface].flux;
        const amrex::IntVect dst      = advection::neighbor_iv(iv, face);

        const bool take_front   = takeFromFront(face);
        const int entry_code    = entryCode(face);

        // Fernando-Clem: only outgoing (negative) fluxes advect from here;
        // incoming mass is the neighbour's export - except at an inflow
        // boundary, where the ghost cell supplies it as a new parcel
        if (net_flux >= 0.0) {
          if ( net_flux > 0.0 && !domain.contains(dst) &&  boundary::is_inflow_boundary_face( dst, domain, is_periodic, lo_bc, hi_bc)) {
            const amrex::IntVect ghost_src = dst;

            ParticleType p_in;
            p_in.id()   = ParticleType::NextID();
            p_in.cpu()  = amrex::ParallelDescriptor::MyProc();

            // Fernando-Clem: primitive incoming state from the bcnormal
            // ghost cell (specific e and Y, unlike the old conservative code)
            const amrex::Real rho_in    = state_arr(ghost_src, rho_indx);
            const amrex::Real one_rho   = 1.0 / rho_in;

            AMREX_D_TERM(p_in.pos(0) = plo[0] + (iv[0] + 0.5) * dx[0];
                         , p_in.pos(1) = plo[1] + (iv[1] + 0.5) * dx[1];
                         , p_in.pos(2) = plo[2] + (iv[2] + 0.5) * dx[2];)

            p_in.rdata(RealData::mass)    = net_flux;
            p_in.rdata(RealData::rho)     = rho_in;
            p_in.rdata(RealData::vol)     = net_flux * one_rho;
            p_in.rdata(RealData::T)       = state_arr(ghost_src, temp_indx);
            p_in.rdata(RealData::eint)    = state_arr(ghost_src, eint_indx) * one_rho;
            p_in.rdata(RealData::xflux)   = net_flux;

            amrex::Real Y_in[NUM_SPECIES];
            for (int sp = 0; sp < NUM_SPECIES; ++sp) {
              Y_in[sp] = state_arr(ghost_src, spec_indx + sp) * one_rho;
              p_in.rdata(RealData::Y0 + sp) = Y_in[sp];
            }
            // Fernando-Clem: incoming element pressure from the EOS of the
            // ghost-cell (bcnormal) state
            amrex::Real P_in = 0.0;
            auto eos = pele::physics::PhysicsType::eos();
            eos.RTY2P(rho_in, p_in.rdata(RealData::T), Y_in, P_in);
            p_in.rdata(RealData::press) = P_in;
            p_in.idata(IntData::lem_index) = 100000;
            // Fernando-Clem: the parcel enters through THIS cell's boundary
            // face: left face -> left end of the line
            if constexpr (
              config::splice_policy == config::SplicePolicy::side_aware) {
              p_in.idata(IntData::entry_side) =
                (face == advection::Face::Left)    ? 1
                : (face == advection::Face::Right) ? 2
                                                   : 3;
            } else {
              p_in.idata(IntData::entry_side) = 1;
            }
            new_particles.push_back(p_in);
            bal.mass_in += net_flux;
          }
          continue;
        }

        const bool exits_through_outflow = boundary::is_outside_domain_on_outflow( dst, domain, is_periodic, lo_bc, hi_bc);

        amrex::Real flux_counter = std::abs(net_flux);

        while (flux_counter > 0.0) {
          if (plist.empty()) {
            amrex::Print() << "CLEM splice: line emptied before the face flux was served\n";
            break;
          }
          const int cand    = take_front ? plist.front() : plist.back();
          ParticleType& p   = particles[cand];

          if (flux_counter > p.rdata(RealData::mass)) {
            // Fernando-Clem: whole parcel leaves through this face
            if (exits_through_outflow) {
              flux_counter      -= p.rdata(RealData::mass);
              bal.mass_out      += p.rdata(RealData::mass);
              p.id()            = -p.id(); // purged at the next Redistribute
            } else {
              const auto wrapped = boundary::wrap_iv(dst, domain, is_periodic);
              const amrex::IntVect& dw = wrapped.first;
              AMREX_D_TERM(p.pos(0) = plo[0] + (dw[0] + 0.5) * dx[0];
                           , p.pos(1) = plo[1] + (dw[1] + 0.5) * dx[1];
                           , p.pos(2) = plo[2] + (dw[2] + 0.5) * dx[2];)
              p.idata(IntData::entry_side)      = entry_code;
              p.idata(IntData::lem_index)       = 100000;
              p.rdata(RealData::xflux)          = std::abs(net_flux);
              flux_counter                      -= p.rdata(RealData::mass);
            }
            // Fernando-Clem: remove the donor from the matching end
            if (take_front) {
              plist.erase(plist.begin());
            } else {
              plist.pop_back();
            }
          } else {
            // Fernando-Clem: split - only flux_counter worth of mass leaves;
            // the remainder keeps its state, volume shrinks with the mass so
            // rho stays intensive
            if (!exits_through_outflow) {
              ParticleType p_out;
              p_out.id()            = ParticleType::NextID();
              p_out.cpu()           = amrex::ParallelDescriptor::MyProc();
              const auto wrapped    = boundary::wrap_iv(dst, domain, is_periodic);
              makeSplitParticle( wrapped.first, plo, dx, flux_counter, net_flux, entry_code, p,p_out);
              new_particles.push_back(p_out);
            } else {
              bal.mass_out += flux_counter;
            }
            p.rdata(RealData::mass)     -= flux_counter;
            p.rdata(RealData::vol)      =  p.rdata(RealData::mass) / p.rdata(RealData::rho);
            flux_counter = 0.0;
          }
        } // while flux_counter
      }   // faces
    }     // iv

    for (const auto& p : new_particles) {
      particle_tile.push_back(p);
    }
  } // pti

  mgr.invalidateCellMap();
  return bal;
}

void
SplicingOperator::reindexAfterSplicing(ClemManager& mgr, const int lev)
{
  BL_PROFILE("clem::SplicingOperator::reindexAfterSplicing()");
  auto& pc = mgr.particleContainer();
  const auto& geom = pc.Geom(lev);
  const auto* dx = geom.CellSize();
  const auto* plo = geom.ProbLo();
  const amrex::Real dx_inner =
    dx[0] / static_cast<amrex::Real>(config::n_lem + 1);

  std::vector<int> plist;

  for (ClemParIter pti(pc, lev); pti.isValid(); ++pti) {
    const int grid_id = pti.index();
    const amrex::Box& bx = pti.validbox();
    auto& particles = pti.GetArrayOfStructs();

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      const auto& cell_list = mgr.particlesInCell(lev, grid_id, iv);
      plist.assign(cell_list.begin(), cell_list.end());

      // Fernando-Clem: Maxwell step 24 splice order - left entries (1) get
      // the lowest indices, residents (0) keep their order in the middle,
      // right/transverse entries (2/3) go to the right end. Within a group
      // the old lem_index ascending preserves the physical line order;
      // stable_sort keeps same-index arrivals deterministic
      std::stable_sort(plist.begin(), plist.end(), [&](int a, int b) {
        auto group = [&](int pidx) {
          const int es = particles[pidx].idata(IntData::entry_side);
          return (es == 1) ? 0 : (es == 0 ? 1 : 2);
        };
        const int ga = group(a);
        const int gb = group(b);
        if (ga != gb) {
          return ga < gb;
        }
        return particles[a].idata(IntData::lem_index) <
               particles[b].idata(IntData::lem_index);
      });

      for (int n = 0; n < static_cast<int>(plist.size()); ++n) {
        ParticleType& p = particles[plist[n]];
        p.idata(IntData::lem_index) = n;
        p.idata(IntData::entry_side) = 0;
        p.rdata(RealData::xflux) = 0.0;
        // Fernando-Clem: reposition only the first n_lem (survivors of the
        // upcoming regrid); moving extras could push them into a neighbour
        // cell and desync them from the map
        if (n < config::n_lem) {
          AMREX_D_TERM(
            p.pos(0) = plo[0] + iv[0] * dx[0] +
                       static_cast<amrex::Real>(n + 1) * dx_inner;
            , p.pos(1) = plo[1] + (iv[1] + 0.5) * dx[1];
            , p.pos(2) = plo[2] + (iv[2] + 0.5) * dx[2];)
        }
      }
    }
  }
}

} // namespace clem

// ================================================================================
// FILE: Source/Clem/ClemDiffusion.H  (original -- edit THERE, not in this bundle)
// ================================================================================
#ifndef CLEMDIFFUSION_H
#define CLEMDIFFUSION_H

#include "PelePhysics.H"

#include "ClemManager.H"
#include "ClemScratch.H"

namespace clem {

// Fernando-Clem: the transport-coefficient block. PeleC owns it (PeleC::
// trans_parms); the host copy is passed in as an argument so the Clem module
// keeps its one-way dependency on PeleC (PeleC includes Clem, never the
// reverse) and stays free of prob.H
using TransParmType = pele::physics::transport::TransParm<
  pele::physics::PhysicsType::eos_type,
  pele::physics::PhysicsType::transport_type>;

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
    // it costs MPI reductions). The line is closed and the element masses are
    // untouched, so both totals are invariant and both numbers are round-off
    // (~1e-15) unless the scheme is broken.
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
    // check for the temperature form: |sum(m h)_after - _before| / sum(m |h|)
    amrex::Real drift_enthalpy = 0.0;
  };

  // Fernando-Clem: advance every line on the level by dt. tparm is the HOST
  // transport-parameter block (the operator runs on the host, over the particle
  // AoS, like the rest of the Clem module)
  static Diagnostics
  diffuse(ClemManager& mgr, int lev, amrex::Real dt, const TransParmType* tparm);

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
};

} // namespace clem
#endif

// ================================================================================
// FILE: Source/Clem/ClemDiffusion.cpp  (original -- edit THERE, not in this bundle)
// ================================================================================
#include <algorithm>
#include <cmath>
#include <limits>

#include <AMReX_ParallelDescriptor.H>

#include "ClemDiffusion.H"
#include "ClemEosUtil.H"

namespace clem {

namespace {

// Fernando-Clem: face value of a transport coefficient. The harmonic mean is
// the series-resistance average: it is the one that does not let a single
// near-zero-conductivity element be bridged by its neighbour, and it degrades
// gracefully to zero flux when either side vanishes
AMREX_FORCE_INLINE amrex::Real
harmonicMean(const amrex::Real a, const amrex::Real b)
{
  const amrex::Real s = a + b;
  return (s > 0.0) ? (2.0 * a * b / s) : 0.0;
}

// Fernando-Clem: below this, a negative Y_k is round-off, not physics. Species
// that are identically zero over half the line (every radical, in the fresh
// gas) round to small negatives constantly; flooring them is right, counting
// them as clipping events is not
constexpr amrex::Real clip_tol = 1.0e-12;

} // namespace

bool
DiffusionOperator::gatherLine(
  LineScratch& lscr,
  const std::vector<int>& plist,
  ClemParticleContainer::ParticleType* particles)
{
  const int n = static_cast<int>(plist.size());
  // Fernando-Clem: a single element has no interior face - nothing to diffuse
  if (n < 2) {
    return false;
  }
  lscr.resize(n);

  // Fernando-Clem: the cell map hands the elements out in AoS order; the
  // operator needs them in LINE order
  lscr.pidx.assign(plist.begin(), plist.begin() + n);
  std::sort(lscr.pidx.begin(), lscr.pidx.begin() + n, [&](const int a, const int b) {
      return particles[a].idata(IntData::lem_index) <  particles[b].idata(IntData::lem_index);
    });

  for (int l = 0; l < n; ++l) {
    const auto& p         = particles[lscr.pidx[l]];
    const amrex::Real m   = p.rdata(RealData::mass);
    const amrex::Real r   = p.rdata(RealData::rho);
    const amrex::Real pr  = p.rdata(RealData::press);
    // Fernando-Clem: dm is a denominator everywhere below, rho scales every
    // coefficient, and press is the constant the isobaric sub-cycle holds - a
    // placeholder element (regrid pads empty cells with them) would poison the
    // whole line, so the line is skipped instead
    if (m <= 0.0 || r <= 0.0 || pr <= 0.0) {
      return false;
    }
    lscr.mass[l]    = m;
    lscr.rho[l]     = r;
    lscr.press[l]   = pr;
    lscr.T[l]       = p.rdata(RealData::T);
    lscr.eint[l]    = p.rdata(RealData::eint);

    for (int k = 0; k < NUM_SPECIES; ++k) {
      lscr.Y[l * NUM_SPECIES + k] = p.rdata(RealData::Y0 + k);
    }
  }
  return true;
}

void
DiffusionOperator::scatterLine(
  const LineScratch& lscr, ClemParticleContainer::ParticleType* particles)
{
  for (int l = 0; l < lscr.n_elem; ++l) {
    auto& p = particles[lscr.pidx[l]];
    p.rdata(RealData::T)    = lscr.T[l];
    p.rdata(RealData::eint) = lscr.eint[l];

    for (int k = 0; k < NUM_SPECIES; ++k) {
      p.rdata(RealData::Y0 + k) = lscr.Y[l * NUM_SPECIES + k];
    }

    // Fernando-Clem: the isobaric sub-cycle held press constant and let rho
    // breathe (rho = rho(press, T, Y)), so press is already the element's true,
    // unchanged value - only density and its volume are written back. vol =
    // mass/rho carries the sub-grid dilatation the LES->SGS coupling needs
    p.rdata(RealData::rho) = lscr.rho[l];
    p.rdata(RealData::vol) = lscr.mass[l] / lscr.rho[l];
  }
}

void
DiffusionOperator::updateProperties(
  LineScratch& lscr, const TransParmType* tparm)
{
  auto eos = pele::physics::PhysicsType::eos();
  auto trans = pele::physics::PhysicsType::transport();

  // Fernando-Clem: only the two coefficients the line equations need
  constexpr bool get_xi = false;
  constexpr bool get_mu = false;
  constexpr bool get_lam = true;
  constexpr bool get_Ddiag = true;
  constexpr bool get_chi = false;

  for (int l = 0; l < lscr.n_elem; ++l) {
    const amrex::Real rho_l = lscr.rho[l];
    const amrex::Real T_l = lscr.T[l];

    amrex::Real Yl[NUM_SPECIES];
    amrex::Real Xl[NUM_SPECIES];
    amrex::Real Cpil[NUM_SPECIES];
    amrex::Real Dl[NUM_SPECIES];
    for (int k = 0; k < NUM_SPECIES; ++k) {
      Yl[k] = lscr.Y[l * NUM_SPECIES + k];
    }

    eos.Y2X(Yl, Xl);
    eos.Y2WBAR(Yl, lscr.wbar[l]);
    eos.RTY2Cp(rho_l, T_l, Yl, lscr.cp[l]); // mixture Cp
    eos.T2Cpi(T_l, Cpil);                   // per-species Cp_k (interspecies term)

    amrex::Real mu        = 0.0;
    amrex::Real xi        = 0.0;
    amrex::Real lam       = 0.0;
    amrex::Real* chi_mix  = nullptr;
    trans.transport(get_xi, get_mu, get_lam, get_Ddiag, get_chi, T_l, rho_l, Yl, Dl, chi_mix,
                    mu, xi, lam, tparm);

    // Fernando-Clem: the mass-coordinate mapping (d/dx = A rho d/dm) applied
    // twice leaves one extra rho on each coefficient; Dl[k] is already
    // Dhat_k = rho D_k W_k/W as returned by PelePhysics
    lscr.rholam[l] = rho_l * lam;

    for (int k = 0; k < NUM_SPECIES; ++k) {
      lscr.X[l * NUM_SPECIES + k]     = Xl[k];
      lscr.cpk[l * NUM_SPECIES + k]   = Cpil[k];
      // Fernando-Clem: rhoD = rho * Dhat_k = rho^2 D_k W_k/W. Simple transport
      // returns Dhat_k = rho D_k W_k/W in Ddiag; the extra rho is the second
      // d/dx -> A rho d/dm mapping. This is the rho*Dhat coefficient of dX_k/dm
      // in CouplingLESwithCLEM.md
      lscr.rhoD[l * NUM_SPECIES + k]  = rho_l * Dl[k];
    }
  }
}

void
DiffusionOperator::computeRates(LineScratch& lscr, const amrex::Real area2)
{
  const int ne = lscr.n_elem;
  const int nf = lscr.n_face; // periodic: nf == ne, face ne-1 wraps ne-1 -> 0

  AMREX_ASSERT(ne == CLEM_NLEM);
  //note that the \Delta m in all particles in CLEM are regridded that is
  // \Delta_m_i = \Delta_m_j for all i,j = 0,...,NUM_CLEM -1;

  // Fernando-Clem: PERIODIC line. Face f sits between element f and element
  // (f+1) mod ne, so it is the "l+1/2" face of element l = f; the last face
  // (f = ne-1) wraps the downstream end back onto element 0 (Maxwell's
  // homogeneous sub-grid domain). Inter-cell transport stays the splicing job
  for (int f = 0; f < nf; ++f) {
    const int lo = f;
    const int hi = (f + 1) % ne;
    // Fernando-Clem: mass-coordinate distance between the two element centres;
    // equals the uniform dm of the doc on a freshly regridded (equal-mass) line
    const amrex::Real odm = 2.0 / (lscr.mass[lo] + lscr.mass[hi]);

    amrex::Real* Sk_f = &lscr.Sk[f * NUM_SPECIES];

    // Fernando-Clem: uncorrected Fickian fluxes F_k = (rho Dhat_k) dX_k/dm and
    // their sum - the sum is exactly the correction (Stefan) velocity that has
    // to be subtracted so the species fluxes close
    amrex::Real Fsum = 0.0;
    for (int k = 0; k < NUM_SPECIES; ++k) {
      const amrex::Real coef  = harmonicMean(lscr.rhoD[lo * NUM_SPECIES + k], lscr.rhoD[hi * NUM_SPECIES + k]);
      const amrex::Real Fk    = coef *(lscr.X[hi * NUM_SPECIES + k] - lscr.X[lo * NUM_SPECIES + k]) * odm;
      Sk_f[k]                 = Fk;
      Fsum                    += Fk;
    }

    // Fernando-Clem: S_k = F_k - Y_k sum_j F_j. Because sum_k Y_k = 1 on both
    // sides of the face, sum_k S_k = 0 identically - that keeps sum_k Y_k = 1 on
    // the line for free. The temperature equation's interspecies transport,
    // sum_k Cp_k S_k dT/dm, is accumulated on the SAME corrected S_k - and it
    // uses the sensible Cp_k, NOT h_k: the formation energy rides with Y and is
    // recovered when e = e(rho, T, Y) is reconstructed after the T update
    const amrex::Real dTf = (lscr.T[hi] - lscr.T[lo]) * odm; // dT/dm at the face
    amrex::Real qint_f = 0.0;
    for (int k = 0; k < NUM_SPECIES; ++k) {
      const amrex::Real Yf   = 0.5 * (lscr.Y[lo * NUM_SPECIES + k] + lscr.Y[hi * NUM_SPECIES + k]);
      Sk_f[k]                -= Yf * Fsum;
      const amrex::Real cpkf = 0.5 * (lscr.cpk[lo * NUM_SPECIES + k] + lscr.cpk[hi * NUM_SPECIES + k]);
      qint_f                 += cpkf * Sk_f[k] * dTf;
    }
    lscr.qint[f] = qint_f;

    // Fernando-Clem: conductive heat flux (rho lambda) dT/dm only - no species-
    // enthalpy flux inside the divergence (that would double-count what qint_f
    // already carries as a gradient product)
    const amrex::Real lam_f = harmonicMean(lscr.rholam[lo], lscr.rholam[hi]);
    lscr.qf[f]              = lam_f * dTf;
  }

  // Fernando-Clem: assemble the rates on the PERIODIC line. Element l has faces
  // l (its l+1/2) and (l-1+ne)%ne (its l-1/2); the wrap makes every element
  // interior. dT/dt = (A^2/Cp)[ d/dm(rho lambda dT/dm) + sum_k Cp_k S_k dT/dm ]:
  // the conduction is a flux divergence (telescopes around the ring), the
  // interspecies term is a pointwise gradient product averaged from the two
  // faces. Species mass telescopes and is conserved; T is non-conservative in
  // energy by construction - the drift is measured in diffuse() (drift_enthalpy)
  for (int l = 0; l < ne; ++l) {
    const int f_hi = l;
    const int f_lo = (l + ne - 1) % ne;
    const amrex::Real cp_l = (lscr.cp[l] > 0.0) ? lscr.cp[l] : 1.0;

    const amrex::Real cond      = (lscr.qf[f_hi] - lscr.qf[f_lo]) / lscr.mass[l];
    const amrex::Real interspec = 0.5 * (lscr.qint[f_hi] + lscr.qint[f_lo]);
    lscr.dTdt[l] = area2 * (cond + interspec) / cp_l;

    const amrex::Real fac = area2 / lscr.mass[l];
    for (int k = 0; k < NUM_SPECIES; ++k) {
      const amrex::Real S_hi          = lscr.Sk[f_hi * NUM_SPECIES + k];
      const amrex::Real S_lo          = lscr.Sk[f_lo * NUM_SPECIES + k];
      lscr.dYdt[l * NUM_SPECIES + k]  = fac * (S_hi - S_lo);
    }
  }
}

amrex::Real
DiffusionOperator::stableSubStep(
  const LineScratch& lscr, const amrex::Real area2)
{
  auto eos    = pele::physics::PhysicsType::eos();
  amrex::Real imw[NUM_SPECIES];
  eos.inv_molecular_weight(imw); // 1/W_k

  // Fernando-Clem: thermal and species transport carry different units, so the
  // two diffusion numbers are formed independently and the tighter one wins
  amrex::Real beta    = 0.0;
  amrex::Real dm_min  = std::numeric_limits<amrex::Real>::max();

  for (int l = 0; l < lscr.n_elem; ++l) {
    dm_min = amrex::min<amrex::Real>(dm_min, lscr.mass[l]);

    // Fernando-Clem: thermal - the T update goes through Cp, hence rho lambda/Cp
    if (lscr.cp[l] > 0.0) {
      beta = amrex::max<amrex::Real>(beta, area2 * lscr.rholam[l] / lscr.cp[l]);
    }

    //Note that  lscr.rhoD = rho * rho D_k W_k/W with \hat{D} = rho D_k W_k/W
    for (int k = 0; k < NUM_SPECIES; ++k) {
      //lets defines Y^{n+1} = G(Y^n) then to be stable, we require an explciit evaluation using Y^n
      //then to do so we defined \rho \hat{D} W/W_k
      const amrex::Real d_eff = lscr.rhoD[l * NUM_SPECIES + k] * lscr.wbar[l] * imw[k];
      beta = amrex::max<amrex::Real>(beta, area2 * d_eff);
    }
  }

  if (beta <= 0.0) {
    return std::numeric_limits<amrex::Real>::max(); // no transport, no limit
  }
  // Fernando-Clem: dm_min (not the mean) keeps the limit conservative on a line
  // left unequal by clem.regrid_int > 1
  return ClemManager::diffusionCfl() * dm_min * dm_min / beta;
}

void
DiffusionOperator::advanceLine(
  LineScratch& lscr, const amrex::Real dt_sub, Diagnostics& diag)
{
  auto eos = pele::physics::PhysicsType::eos();

  for (int l = 0; l < lscr.n_elem; ++l) {
    // Fernando-Clem: T is the transported (primary) thermal variable
    lscr.T[l] += dt_sub * lscr.dTdt[l];

    amrex::Real Yl[NUM_SPECIES];
    amrex::Real Ysum = 0.0;
    amrex::Real worst_neg = 0.0;
    for (int k = 0; k < NUM_SPECIES; ++k) {
      amrex::Real Yk = lscr.Y[l * NUM_SPECIES + k] + dt_sub * lscr.dYdt[l * NUM_SPECIES + k];
      if (Yk < 0.0) {
        worst_neg = amrex::max<amrex::Real>(worst_neg, -Yk);
        Yk = 0.0;
      }
      Yl[k] = Yk;
      Ysum += Yk;
    }
    
    //This is only if negative values exist
    if (worst_neg > 0.0) {
      // Fernando-Clem: sum_k Y_k = 1 is analytic here (sum_k S_k = 0), so the
      // renormalisation only ever undoes what the flooring just added
      const amrex::Real onorm = (Ysum > 0.0) ? 1.0 / Ysum : 0.0;
      for (int k = 0; k < NUM_SPECIES; ++k) {
        Yl[k] *= onorm;
      }
      diag.max_neg_Y = amrex::max<amrex::Real>(diag.max_neg_Y, worst_neg);
      if (worst_neg > clip_tol) {
        ++diag.n_clipped;
      }
    }

    for (int k = 0; k < NUM_SPECIES; ++k) {
      lscr.Y[l * NUM_SPECIES + k] = Yl[k];
    }

    // Fernando-Clem: isobaric closure - density and internal energy follow from
    // the new (T, Y) at the HELD element pressure (PYT2RE returns both). The
    // line breathes (rho changes) while the mass grid stays fixed
    amrex::Real rho_new = 0.0;
    amrex::Real e_new   = 0.0;
    eos.PYT2RE(lscr.press[l], Yl, lscr.T[l], rho_new, e_new);
    lscr.rho[l]  = rho_new;
    lscr.eint[l] = e_new;
  }
}

DiffusionOperator::Diagnostics
DiffusionOperator::diffuse(
  ClemManager& mgr,
  const int lev,
  const amrex::Real dt,
  const TransParmType* tparm)
{
  BL_PROFILE("clem::DiffusionOperator::diffuse()");
  AMREX_ALWAYS_ASSERT(mgr.isDefined());
  AMREX_ALWAYS_ASSERT(tparm != nullptr);

  auto& pc = mgr.particleContainer();
  if (!mgr.IsCellMapValid()) {
    mgr.rebuildCellMap(lev);
  }

  // Fernando-Clem: the LEM line is x-aligned and spans the cell, so its
  // cross-section is the cell volume over the line length, A = V_cell/dx.
  // (3D: dy*dz; 2D: dy; 1D: 1 - i.e. exactly AMReX's unit-depth convention, the
  // same one writeBackToLes uses to form the cell volume)
  const auto* dx = pc.Geom(lev).CellSize();
  amrex::Real cell_vol = 1.0;
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    cell_vol *= dx[d];
  }
  const amrex::Real area    = cell_vol / dx[0];
  const amrex::Real area2   = area * area;

  const int max_sub = ClemManager::diffusionMaxSubsteps();

  Diagnostics diag;
  // Fernando-Clem: ONE arena for the whole level - gathered, worked on and
  // scattered back cell by cell (see ClemScratch.H)
  auto& lscr = LineScratch::get();

  // Fernando-Clem: enthalpy-drift accumulators. Constant-pressure diffusion
  // conserves H = sum(m h), not internal energy; the T-transport is non-
  // conservative, so H_after - H_before measures the scheme's energy drift
  amrex::Real H_before = 0.0;
  amrex::Real H_after  = 0.0;
  amrex::Real H_abs    = 0.0;

  for (ClemParIter pti(pc, lev); pti.isValid(); ++pti) {
    const int grid_id       = pti.index();
    const amrex::Box& bx    = pti.validbox();
    auto* particles         = pti.GetArrayOfStructs().data();

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      const auto& plist = mgr.particlesInCell(lev, grid_id, iv);
      if (!gatherLine(lscr, plist, particles)) {
        ++diag.n_skipped;
        continue;
      }

      // Fernando-Clem: the properties of the first sub-step double as the ones
      // the stability limit is measured on, so they are computed once here and
      // the sub-step loop refreshes them only from its second pass on
      updateProperties(lscr, tparm);

      // Fernando-Clem: line enthalpy before diffusion, h = e + p/rho
      for (int l = 0; l < lscr.n_elem; ++l) {
        H_before += lscr.mass[l] * (lscr.eint[l] + lscr.press[l] / lscr.rho[l]);
      }

      const amrex::Real dt_cfl = stableSubStep(lscr, area2);
      int n_sub = 1;
      if (dt_cfl < dt) {
        const amrex::Real want = std::ceil(dt / dt_cfl);
        if (want >= static_cast<amrex::Real>(max_sub)) {
          // Fernando-Clem: the cap wins, so this line is about to be integrated
          // ABOVE its stability limit. Silently capping is how an unstable line
          // turns into a NaN three steps later, so it is counted and shouted
          n_sub = max_sub;
          ++diag.n_capped;
        } else {
          n_sub = static_cast<int>(want);
        }
      }
      const amrex::Real dt_sub = dt / static_cast<amrex::Real>(n_sub);
      diag.max_substeps = amrex::max<int>(diag.max_substeps, n_sub);
      // Fernando-Clem: accumulate the TOTAL work (summed across ranks below):
      // one unit = one full-line computeRates+advanceLine sweep
      diag.total_substeps += static_cast<amrex::Long>(n_sub);
      ++diag.n_lines;

      for (int s = 0; s < n_sub; ++s) {
        if (s > 0) {
          // Fernando-Clem: coefficients AND driving forces are re-derived from
          // the current (T, Y) so they stay consistent as the line evolves
          updateProperties(lscr, tparm);
        }
        computeRates(lscr, area2);
        advanceLine(lscr, dt_sub, diag);
      }

      // Fernando-Clem: line enthalpy after diffusion (rho has breathed)
      for (int l = 0; l < lscr.n_elem; ++l) {
        const amrex::Real h = lscr.eint[l] + lscr.press[l] / lscr.rho[l];
        H_after += lscr.mass[l] * h;
        H_abs   += lscr.mass[l] * std::abs(h);
      }

      scatterLine(lscr, particles);
    } // iv
  }   // pti

  amrex::ParallelDescriptor::ReduceIntMax(diag.max_substeps);
  amrex::ParallelDescriptor::ReduceLongSum(diag.total_substeps);
  amrex::ParallelDescriptor::ReduceLongSum(diag.n_lines);
  amrex::ParallelDescriptor::ReduceLongSum(diag.n_skipped);
  amrex::ParallelDescriptor::ReduceLongSum(diag.n_clipped);
  amrex::ParallelDescriptor::ReduceLongSum(diag.n_capped);
  amrex::ParallelDescriptor::ReduceRealMax(diag.max_neg_Y);

  // Fernando-Clem: enthalpy drift over the whole level (the constant-p invariant
  // the T-transport should preserve); non-zero is the scheme's energy error
  amrex::ParallelDescriptor::ReduceRealSum(H_before);
  amrex::ParallelDescriptor::ReduceRealSum(H_after);
  amrex::ParallelDescriptor::ReduceRealSum(H_abs);
  diag.drift_enthalpy = (H_abs > 0.0) ? std::abs(H_after - H_before) / H_abs : 0.0;

  // Fernando-Clem: the count is already reduced, so every rank holds it - only
  // one of them should say so
  if (diag.n_capped > 0 && amrex::ParallelDescriptor::IOProcessor()) {
    amrex::Warning(
      "CLEM diffusion: clem.diffusion_max_substeps was hit - those lines were "
      "integrated above their stability limit and WILL go negative and blow "
      "up. Raise clem.diffusion_max_substeps (cost grows linearly) or coarsen "
      "the line; do not ignore this.");
  }
  return diag;
}

} // namespace clem

// ================================================================================
// FILE: Source/Clem/ClemLesCoupling.H  (original -- edit THERE, not in this bundle)
// ================================================================================
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

// ================================================================================
// FILE: Source/Clem/ClemLesCoupling.cpp  (original -- edit THERE, not in this bundle)
// ================================================================================
#include <cmath>

#include <AMReX_ParallelDescriptor.H>

#include "PelePhysics.H"

#include "ClemLesCoupling.H"

namespace clem {

void
LesCoupling::isentropicPressureUpdate(
  ClemManager& mgr,
  const int lev,
  const amrex::MultiFab& s_new,
  const int rho_indx,
  const int temp_indx,
  const int spec_indx)
{
  BL_PROFILE("clem::LesCoupling::isentropicPressureUpdate()");
  auto& pc = mgr.particleContainer();
  AMREX_ALWAYS_ASSERT(s_new.boxArray() == pc.ParticleBoxArray(lev));
  if (!mgr.IsCellMapValid()) {
    mgr.rebuildCellMap(lev);
  }
  auto eos = pele::physics::PhysicsType::eos();

  for (ClemParIter pti(pc, lev); pti.isValid(); ++pti) {
    const int grid_id    = pti.index();
    const amrex::Box& bx = pti.validbox();
    const auto& sn_arr   = s_new.const_array(grid_id);
    auto& particles      = pti.GetArrayOfStructs();

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      const auto& plist = mgr.particlesInCell(lev, grid_id, iv);
      if (plist.empty()) {
        continue;
      }

      // Fernando-Clem: target pressure = resolved cell pressure of the new
      // LES state, computed once per cell
      amrex::Real P_new = 0.0;
      {
        amrex::Real Y_les[NUM_SPECIES];
        const amrex::Real r_n   = sn_arr(iv, rho_indx);
        const amrex::Real onr_n = 1.0 / r_n;
        for (int sp = 0; sp < NUM_SPECIES; ++sp) {
          Y_les[sp] = sn_arr(iv, spec_indx + sp) * onr_n;
        }
        eos.RTY2P(r_n, sn_arr(iv, temp_indx), Y_les, P_new);
      }

      for (const int pidx : plist) {
        auto& p                   = particles[pidx];
        const amrex::Real rho_old = p.rdata(RealData::rho);
        const amrex::Real T_old   = p.rdata(RealData::T);
        // Fernando-Clem: reference pressure = the element's OWN pressure
        // (Maxwell element state p_k)
        const amrex::Real P_old   = p.rdata(RealData::press);
        if (rho_old <= 0.0 || P_old <= 0.0) {
          continue; // placeholder element
        }

        // Fernando-Clem: mass fractions are stored primitive and are
        // invariant under the isentropic compression
        amrex::Real Y[NUM_SPECIES];
        for (int sp = 0; sp < NUM_SPECIES; ++sp) {
          Y[sp] = p.rdata(RealData::Y0 + sp);
        }

        // Fernando-Clem: local isentropic exponent so the compression
        // tracks the real mixture instead of a hard-coded gamma
        amrex::Real gamma = 0.0;
        eos.RTY2G(rho_old, T_old, Y, gamma);

        const amrex::Real Pratio  = P_new / P_old;
        const amrex::Real rho_fac = std::pow(Pratio, 1.0 / gamma);
        const amrex::Real T_new   = T_old * std::pow(Pratio, (gamma - 1.0) / gamma);
        const amrex::Real rho_new = rho_old * rho_fac;

        p.rdata(RealData::T)     = T_new;
        p.rdata(RealData::rho)   = rho_new;
        p.rdata(RealData::press) = P_new;
        // Fernando-Clem: mass fixed under compression -> volume shrinks
        p.rdata(RealData::vol) /= rho_fac;

        // Fernando-Clem: DIFFERENCE-form energy update,
        //   e += e(T_new, Y) - e(T_old, Y),
        // so any constant T-representation bias cancels and the stored
        // energy changes only by the isentropic work (no per-step ratchet)
        amrex::Real e_Told = 0.0;
        amrex::Real e_Tnew = 0.0;
        eos.RTY2E(rho_old, T_old, Y, e_Told);
        eos.RTY2E(rho_new, T_new, Y, e_Tnew);
        p.rdata(RealData::eint) += e_Tnew - e_Told;
      }
    }
  }
}

amrex::Real
LesCoupling::renormalizeEnsembleVolume(ClemManager& mgr, const int lev)
{
  BL_PROFILE("clem::LesCoupling::renormalizeEnsembleVolume()");
  auto& pc = mgr.particleContainer();
  if (!mgr.IsCellMapValid()) {
    mgr.rebuildCellMap(lev);
  }
  auto eos = pele::physics::PhysicsType::eos();

  const auto* dx = pc.Geom(lev).CellSize();
  amrex::Real cell_vol = 1.0;
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    cell_vol *= dx[d];
  }

  amrex::Real max_defect = 0.0;

  for (ClemParIter pti(pc, lev); pti.isValid(); ++pti) {
    const int grid_id    = pti.index();
    const amrex::Box& bx = pti.validbox();
    auto& particles      = pti.GetArrayOfStructs();

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      const auto& plist = mgr.particlesInCell(lev, grid_id, iv);
      if (plist.empty()) {
        continue;
      }

      // Fernando-Clem: current line volume (placeholders carry vol = 0)
      amrex::Real v_line = 0.0;
      for (const int pidx : plist) {
        v_line += particles[pidx].rdata(RealData::vol);
      }
      if (v_line <= 0.0) {
        continue;
      }

      // Fernando-Clem: uniform volume ratio putting the line back into the cell
      const amrex::Real f = cell_vol / v_line;
      max_defect = amrex::max<amrex::Real>(max_defect, std::abs(f - 1.0));
      if (f == 1.0) {
        continue;
      }

      for (const int pidx : plist) {
        auto& p                   = particles[pidx];
        const amrex::Real rho_old = p.rdata(RealData::rho);
        const amrex::Real T_old   = p.rdata(RealData::T);
        const amrex::Real P_old   = p.rdata(RealData::press);
        if (rho_old <= 0.0 || P_old <= 0.0) {
          continue; // placeholder element
        }

        amrex::Real Y[NUM_SPECIES];
        for (int sp = 0; sp < NUM_SPECIES; ++sp) {
          Y[sp] = p.rdata(RealData::Y0 + sp);
        }
        amrex::Real gamma = 0.0;
        eos.RTY2G(rho_old, T_old, Y, gamma);

        // Fernando-Clem: isentropic rescale by the SAME f for every element:
        // vol -> f vol, rho -> rho/f, T V^(g-1) = const, p V^g = const
        const amrex::Real rho_new = rho_old / f;
        const amrex::Real T_new   = T_old * std::pow(f, -(gamma - 1.0));
        const amrex::Real P_new   = P_old * std::pow(f, -gamma);

        p.rdata(RealData::rho)   = rho_new;
        p.rdata(RealData::T)     = T_new;
        p.rdata(RealData::press) = P_new;
        p.rdata(RealData::vol)   = p.rdata(RealData::mass) / rho_new;

        // Fernando-Clem: DIFFERENCE-form energy update (same convention as the
        // isentropic pressure update) - only the compression work enters
        amrex::Real e_old = 0.0;
        amrex::Real e_new = 0.0;
        eos.RTY2E(rho_old, T_old, Y, e_old);
        eos.RTY2E(rho_new, T_new, Y, e_new);
        //
        p.rdata(RealData::eint) += e_new - e_old;
      }
    }
  }

  amrex::ParallelDescriptor::ReduceRealMax(max_defect);
  return max_defect;
}

void
LesCoupling::snapshotFilteredMean(
  ClemManager& mgr,
  const int lev,
  amrex::MultiFab& out,
  const bool conserved)
{
  BL_PROFILE("clem::LesCoupling::snapshotFilteredMean()");
  auto& pc = mgr.particleContainer();
  AMREX_ALWAYS_ASSERT(out.boxArray() == pc.ParticleBoxArray(lev));
  AMREX_ALWAYS_ASSERT(out.nComp() == NUM_SPECIES + 1);
  if (!mgr.IsCellMapValid()) {
    mgr.rebuildCellMap(lev);
  }

  // Fernando-Clem: conserved form divides the mass-weighted sums by the cell
  // volume (giving rho Y_k, rho e); ratio form divides by the cell mass sum
  const auto* dx = pc.Geom(lev).CellSize();
  amrex::Real cell_vol = 1.0;
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    cell_vol *= dx[d];
  }
  const amrex::Real one_vol = 1.0 / cell_vol;

  // Fernando-Clem: empty cells stay 0 so pre and post cancel exactly there
  out.setVal(0.0);

  for (ClemParIter pti(pc, lev); pti.isValid(); ++pti) {
    const int grid_id     = pti.index();
    const amrex::Box& bx  = pti.validbox();
    const auto& o_arr     = out.array(grid_id);
    auto& particles       = pti.GetArrayOfStructs();

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      const auto& plist = mgr.particlesInCell(lev, grid_id, iv);
      if (plist.empty()) {
        continue;
      }

      // Fernando-Clem: filtered (mass-weighted) mean of the cell ensemble
      amrex::Real m_sum               = 0.0;
      amrex::Real me_sum              = 0.0;
      amrex::Real mY_sum[NUM_SPECIES] = {0.0};
      for (const int pidx : plist) {
        const auto& p       = particles[pidx];
        const amrex::Real m = p.rdata(RealData::mass);
        m_sum  += m;
        me_sum += m * p.rdata(RealData::eint);
        for (int sp = 0; sp < NUM_SPECIES; ++sp) {
          mY_sum[sp] += m * p.rdata(RealData::Y0 + sp);
        }
      }
      if (m_sum <= 0.0) {
        continue;
      }
      const amrex::Real scale = conserved ? one_vol : (1.0 / m_sum);
      for (int sp = 0; sp < NUM_SPECIES; ++sp) {
        o_arr(iv, sp) = mY_sum[sp] * scale;
      }
      o_arr(iv, NUM_SPECIES) = me_sum * scale;
    }
  }
}

amrex::Real
LesCoupling::writeBackProjected(
  ClemManager& mgr,
  const int lev,
  amrex::MultiFab& s_new,
  const amrex::MultiFab& s_old,
  const amrex::MultiFab& cons_pre,
  const amrex::MultiFab& cons_post,
  const int rho_indx,
  const int mom_indx,
  const int eden_indx,
  const int temp_indx,
  const int eint_indx,
  const int spec_indx)
{
  BL_PROFILE("clem::LesCoupling::writeBackProjected()");
  auto eos = pele::physics::PhysicsType::eos();
  auto& pc = mgr.particleContainer();
  AMREX_ALWAYS_ASSERT(s_new.boxArray() == pc.ParticleBoxArray(lev));
  AMREX_ALWAYS_ASSERT(cons_pre.nComp() == NUM_SPECIES + 1);
  if (!mgr.IsCellMapValid()) {
    mgr.rebuildCellMap(lev);
  }

  amrex::Real max_proj = 0.0; // max |rho_new/sum - 1| = the realizability residual

  for (ClemParIter pti(pc, lev); pti.isValid(); ++pti) {
    const int grid_id    = pti.index();
    const amrex::Box& bx = pti.validbox();
    const auto& sn       = s_new.array(grid_id);
    const auto& so       = s_old.const_array(grid_id);
    const auto& cpre     = cons_pre.const_array(grid_id);
    const auto& cpost    = cons_post.const_array(grid_id);

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      const auto& plist = mgr.particlesInCell(lev, grid_id, iv);
      if (plist.empty()) {
        continue;
      }
      const amrex::Real rho_new = sn(iv, rho_indx);
      if (rho_new <= 0.0) {
        continue;
      }

      // Fernando-Clem: (rho Y_k)^n_LES + CLEM conserved delta (splice transport
      // + diffusion + reaction), floored to the realizable simplex
      amrex::Real rhoY[NUM_SPECIES];
      amrex::Real sum = 0.0;
      for (int sp = 0; sp < NUM_SPECIES; ++sp) {
        amrex::Real v = so(iv, spec_indx + sp) + (cpost(iv, sp) - cpre(iv, sp));
        if (v < 0.0) {
          v = 0.0;
        }
        rhoY[sp] = v;
        sum += v;
      }
      if (sum <= 0.0) {
        continue;
      }

      // Fernando-Clem: PROJECT onto sum_k rho Y_k = rho_LES^{n+1} (mass
      // realizability); the rescale factor differs from 1 only by the splicing
      // mass-quantization residual (~rho_err)
      const amrex::Real fac = rho_new / sum;
      max_proj = amrex::max<amrex::Real>(max_proj, std::abs(fac - 1.0));

      for (int sp = 0; sp < NUM_SPECIES; ++sp) {
        sn(iv, spec_indx + sp) = rhoY[sp] * fac;
      }

      // Fernando-Clem: energy treatment (clem.couple_back_energy mode).
      const int e_mode = ClemManager::coupleBackEnergyMode();
      if (e_mode == 1) {
        // mode 1 - conserved delta, the SAME scheme as the species:
        // (rho e)^{n+1} = (rho e)^n_LES + Delta_CLEM. cons_* are CONSERVED
        // snapshots, so component NUM_SPECIES is already (rho e)_filtered =
        // sum(m_p e_p)/V_cell [erg/cm^3]; the delta is added DIRECTLY (no rho
        // factor). The CLEM delta carries the isentropic pressure work (mirror
        // of the LES -p div u) + the spliced energy transport (mirror of the
        // convective flux). Measured NEUTRAL in pure advection.
        amrex::Real old_ke = 0.0;
        for (int m = 0; m < 3; ++m) {
          old_ke += so(iv, mom_indx + m) * so(iv, mom_indx + m);
        }
        old_ke *= 0.5 / so(iv, rho_indx);

        const amrex::Real rho_e_old = so(iv, eden_indx) - old_ke;
        const amrex::Real rho_e = rho_e_old + (cpost(iv, NUM_SPECIES) - cpre(iv, NUM_SPECIES));

        sn(iv, eint_indx) = rho_e;
        amrex::Real new_ke = 0.0;
        for (int m = 0; m < 3; ++m) {
          new_ke += sn(iv, mom_indx + m) * sn(iv, mom_indx + m);
        }
        new_ke *= 0.5 / rho_new;
        sn(iv, eden_indx) = rho_e + new_ke;
      } else if (e_mode == 2) {
        // mode 2 - keep-T reconciliation: the LES temperature (temp_indx, made
        // EOS-consistent by computeTemp BEFORE the CLEM stage) is preserved
        // under the composition change: rho e := rho * e(T_pre, Y_new). The
        // subsequent computeTemp then recovers T_pre exactly, so the conduction
        // operator never sees the EOS T-kick of the write-back dY (the
        // velocity-deficit mechanism). The injected energy is the formation-
        // energy of dY at fixed T - non-conservative; monitor for pumping.
        const amrex::Real T_pre = sn(iv, temp_indx);
        amrex::Real Ynew[NUM_SPECIES];
        const amrex::Real onr = 1.0 / rho_new;
        for (int sp = 0; sp < NUM_SPECIES; ++sp) {
          Ynew[sp] = sn(iv, spec_indx + sp) * onr; // post-projection composition
        }
        amrex::Real e_keepT = 0.0;
        eos.RTY2E(rho_new, T_pre, Ynew, e_keepT);
        const amrex::Real rho_e = rho_new * e_keepT;

        sn(iv, eint_indx) = rho_e;
        amrex::Real new_ke = 0.0;
        for (int m = 0; m < 3; ++m) {
          new_ke += sn(iv, mom_indx + m) * sn(iv, mom_indx + m);
        }
        new_ke *= 0.5 / rho_new;
        sn(iv, eden_indx) = rho_e + new_ke;
      }


    }
  }

  amrex::ParallelDescriptor::ReduceRealMax(max_proj);
  return max_proj;
}

amrex::Real
LesCoupling::writeBackToLes(
  ClemManager& mgr,
  const int lev,
  amrex::MultiFab& s_new,
  amrex::MultiFab& s_old,
  const amrex::MultiFab& sgs_incr,
  const int rho_indx,
  const int mom_indx,
  const int eden_indx,
  const int temp_indx,
  const int eint_indx,
  const int spec_indx)
{
  BL_PROFILE("clem::LesCoupling::writeBackToLes()");
  amrex::ignore_unused(temp_indx);
  auto& pc = mgr.particleContainer();
  AMREX_ALWAYS_ASSERT(s_new.boxArray() == pc.ParticleBoxArray(lev));
  AMREX_ALWAYS_ASSERT(sgs_incr.nComp() == NUM_SPECIES + 1);
  if (!mgr.IsCellMapValid()) {
    mgr.rebuildCellMap(lev);
  }

  const auto* dx = pc.Geom(lev).CellSize();
  amrex::Real cell_vol = 1.0;
  
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    cell_vol *= dx[d];
  }
  const amrex::Real one_vol = 1.0 / cell_vol;

  amrex::Real rho_err = 0.0;

  for (ClemParIter pti(pc, lev); pti.isValid(); ++pti) {
    const int grid_id     = pti.index();
    const amrex::Box& bx  = pti.validbox();
    const auto& sn_arr    = s_new.array(grid_id);
    const auto& so_arr    = s_old.array(grid_id);
    const auto& inc_arr   = sgs_incr.const_array(grid_id);
    auto& particles       = pti.GetArrayOfStructs();

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      const auto& plist = mgr.particlesInCell(lev, grid_id, iv);
      if (plist.empty()) {
        continue;
      }

      // Fernando-Clem: filtered mass for the density coupling diagnostic
      amrex::Real m_sum = 0.0;
      for (const int pidx : plist) {
        m_sum += particles[pidx].rdata(RealData::mass);
      }
      if (m_sum <= 0.0) {
        continue;
      }

      // Fernando-Clem: coupling quality - the filtered density must match
      // the LES density up to splicing round-off (same mass fluxes)
      const amrex::Real rho_les = sn_arr(iv, rho_indx);
      const amrex::Real rho_bar = m_sum * one_vol;
      rho_err = amrex::max<amrex::Real>(rho_err, std::abs(rho_bar - rho_les) / rho_les);

      // Fernando-Clem: INCREMENT feedback. Add only what the subgrid processes
      // (diffusion + reaction) changed the filtered composition by; the LES
      // keeps its own 2nd-order resolved transport. sum_k of the increment is
      // zero (mass conserved by the processes), so sum_k rhoY_k stays = rho_LES
      // and no spurious pressure/velocity bulge is injected. In pure advection
      // the increment is identically zero -> the LES is left untouched.
      for (int sp = 0; sp < NUM_SPECIES; ++sp) {
        sn_arr(iv, spec_indx + sp) += rho_les * inc_arr(iv, sp);
      }

      // Fernando-Clem: energy increment (optional). The isentropic LES->SGS
      // pressure work is NOT in this bracket (pre is taken after it), so only
      // the process (reaction/diffusion) energy change is fed back - the same
      // conservative source-only idea as React.cpp. rho*e -> rho*e + rho_LES*de
      //The SGS model only add energy due to reaction so this cannnot be added.
      if (ClemManager::coupleBackEnergy()) {
        const amrex::Real xvel_old = so_arr(iv,mom_indx ) / sn_arr(iv, rho_indx);
        const amrex::Real yvel_old = so_arr(iv,mom_indx ) / sn_arr(iv, rho_indx);
        const amrex::Real zvel_old = so_arr(iv,mom_indx ) / sn_arr(iv, rho_indx);
        const amrex::Real old_ke = sn_arr(iv, rho_indx) * (xvel_old * xvel_old +  yvel_old * yvel_old + zvel_old * zvel_old) * 0.5;


        const amrex::Real rho_e_old = sn_arr(iv, eden_indx) - old_ke;
        //const amrex::Real rho_e = sn_arr(iv, eint_indx) + rho_les * inc_arr(iv, NUM_SPECIES);
        const amrex::Real rho_e = rho_e_old +  rho_les * inc_arr(iv, NUM_SPECIES);

        sn_arr(iv, eint_indx) = rho_e;
        amrex::Real new_ke = 0.0;
        for (int m = 0; m < 3; ++m) {
          new_ke += sn_arr(iv, mom_indx + m) * sn_arr(iv, mom_indx + m);
        }
        new_ke *= 0.5 / rho_les;
        sn_arr(iv, eden_indx) = rho_e + new_ke;
      }
    }
  }

  amrex::ParallelDescriptor::ReduceRealMax(rho_err);
  return rho_err;
}

amrex::Real
LesCoupling::conservedConsistency(
  ClemManager& mgr,
  const int lev,
  const amrex::MultiFab& cons_pre,
  const amrex::MultiFab& cons_post,
  const amrex::MultiFab& s_new,
  const amrex::MultiFab& s_old,
  const int rho_indx,
  const int eint_indx,
  const int spec_indx)
{
  BL_PROFILE("clem::LesCoupling::conservedConsistency()");
  amrex::ignore_unused(rho_indx);
  // Fernando-Clem: LES-vs-SGS conserved instrument (see header). Per cell the
  // two machineries act on the same conserved filtered fields; the per-step gap
  //   d(rho phi)_CLEM - d(rho phi)_LES,  phi in {Y_k, e}
  // localizes the CLEM bookkeeping error by CHANNEL (species mass vs energy),
  // and the pre-offset |cons_pre - s_old| tracks the accumulated drift.
  auto& pc = mgr.particleContainer();
  if (!mgr.IsCellMapValid()) {
    mgr.rebuildCellMap(lev);
  }
  const auto* plo = pc.Geom(lev).ProbLo();
  const auto* dx  = pc.Geom(lev).CellSize();

  // species channel [g/cm^3] and energy channel [erg/cm^3] tracked separately
  amrex::Real gapY = 0.0, dY_les_max = 0.0, offY = 0.0, xY = -1.0;
  amrex::Real gapE = 0.0, dE_les_max = 0.0, offE = 0.0, xE = -1.0;

  for (ClemParIter pti(pc, lev); pti.isValid(); ++pti) {
    const int grid_id    = pti.index();
    const amrex::Box& bx = pti.validbox();
    const auto& pre      = cons_pre.const_array(grid_id);
    const auto& post     = cons_post.const_array(grid_id);
    const auto& sn_arr   = s_new.const_array(grid_id);
    const auto& so_arr   = s_old.const_array(grid_id);

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      const auto& plist = mgr.particlesInCell(lev, grid_id, iv);
      if (plist.empty()) {
        continue;
      }

      // species mass channel
      for (int sp = 0; sp < NUM_SPECIES; ++sp) {
        const amrex::Real d_les  = sn_arr(iv, spec_indx + sp) - so_arr(iv, spec_indx + sp);
        const amrex::Real d_clem = post(iv, sp) - pre(iv, sp);
        const amrex::Real gap    = std::abs(d_clem - d_les);
        dY_les_max = amrex::max<amrex::Real>(dY_les_max, std::abs(d_les));
        offY       = amrex::max<amrex::Real>(offY, std::abs(pre(iv, sp) - so_arr(iv, spec_indx + sp)));
        if (gap > gapY) {
          gapY = gap;
          xY   = plo[0] + (iv[0] + 0.5) * dx[0];
        }
      }

      // energy channel (rho e; eint_indx carries the LES internal energy)
      {
        const amrex::Real d_les  = sn_arr(iv, eint_indx) - so_arr(iv, eint_indx);
        const amrex::Real d_clem = post(iv, NUM_SPECIES) - pre(iv, NUM_SPECIES);
        const amrex::Real gap    = std::abs(d_clem - d_les);
        dE_les_max = amrex::max<amrex::Real>(dE_les_max, std::abs(d_les));
        offE       = amrex::max<amrex::Real>(offE, std::abs(pre(iv, NUM_SPECIES) - so_arr(iv, eint_indx)));
        if (gap > gapE) {
          gapE = gap;
          xE   = plo[0] + (iv[0] + 0.5) * dx[0];
        }
      }
    }
  }

  // Fernando-Clem: max-loc reductions (x >= 0 in this domain masks the losers)
  const amrex::Real myY = gapY, myE = gapE;
  amrex::ParallelDescriptor::ReduceRealMax(gapY);
  amrex::ParallelDescriptor::ReduceRealMax(gapE);
  amrex::ParallelDescriptor::ReduceRealMax(dY_les_max);
  amrex::ParallelDescriptor::ReduceRealMax(dE_les_max);
  amrex::ParallelDescriptor::ReduceRealMax(offY);
  amrex::ParallelDescriptor::ReduceRealMax(offE);
  if (myY < gapY) { xY = -1.0; }
  if (myE < gapE) { xE = -1.0; }
  amrex::ParallelDescriptor::ReduceRealMax(xY);
  amrex::ParallelDescriptor::ReduceRealMax(xE);

  if (ClemManager::verbose() > 0) {
    const amrex::Real relY = (dY_les_max > 0.0) ? gapY / dY_les_max : 0.0;
    const amrex::Real relE = (dE_les_max > 0.0) ? gapE / dE_les_max : 0.0;
    amrex::Print() << "CLEM LESvsSGS[rhoY]: gap = " << gapY
                   << " (rel " << relY << ") at x = " << xY
                   << ", |dLES| = " << dY_les_max
                   << ", sync-off = " << offY << '\n'
                   << "CLEM LESvsSGS[rhoE]: gap = " << gapE
                   << " (rel " << relE << ") at x = " << xE
                   << ", |dLES| = " << dE_les_max
                   << ", sync-off = " << offE << '\n';
  }
  return gapY;
}

} // namespace clem

// ================================================================================
// FILE: Source/Clem/ClemAlgorithm.H  (original -- edit THERE, not in this bundle)
// ================================================================================
#ifndef CLEMALGORITHM_H
#define CLEMALGORITHM_H

#include "ClemDiffusion.H"
#include "ClemManager.H"
#include "ClemSplicing.H"

namespace clem {

// Fernando-Clem: per-step subgrid algorithm (alpaca-style module): the stage
// sequence of one CLEM step is assembled here from the operator modules, not
// inside PeleC. Current composition (Maxwell chain):
//   1. scale the accumulated face mass fluxes (time centering)
//   2. LesCoupling::isentropicPressureUpdate - LES->SGS pressure work
//   3. Regridder - restore n_lem equal-mass elements per cell
//   4. DiffusionOperator - molecular diffusion along the line (clem.do_diffusion)
//   5. (stirring -> reaction: inserted here when implemented)
//   6. SplicingOperator - inter-cell transport with the LES mass fluxes
//   7. LesCoupling::writeBackToLes - SGS->LES filtered (rho Y_k, rho e)
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
  };

  // Fernando-Clem: advance the subgrid by dt. flux_scale converts the
  // accumulated MOL-stage fluxes to the time-centered flux (0.5 for the
  // two-stage Heun scheme); state_bc must be FillPatched with >= 1 ghost
  // cell (bcnormal inflow state); s_new is the LES state at t^{n+1} (its
  // pressure drives the isentropic stage; its composition/energy are
  // overwritten by the write-back when clem.couple_back is on); nstep drives
  // the periodic particle dumps; tparm is PeleC's HOST transport-parameter
  // block, needed by the diffusion stage
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
    const TransParmType* tparm);
};

} // namespace clem
#endif

// ================================================================================
// FILE: Source/Clem/ClemAlgorithm.cpp  (original -- edit THERE, not in this bundle)
// ================================================================================
#include <cmath>

#include <AMReX_Utility.H>

#include "ClemAlgorithm.H"
#include "ClemLesCoupling.H"
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
  // energy among the elements of a cell and leaves the element masses alone
  if (ClemManager::doDiffusion()) {
    const bool check_conservation = (ClemManager::verbose() > 1);
    ConservedTotals before;
    if (check_conservation) {
      before = conservedTotals(mgr, lev);
    }

    diag.diffusion = DiffusionOperator::diffuse(mgr, lev, dt, tparm);

    if (check_conservation) {
      const ConservedTotals after = conservedTotals(mgr, lev);
      diag.diffusion.drift_energy = driftEnergy(before, after);
      diag.diffusion.drift_species = driftSpecies(before, after);
    }
  }

  // Fernando-Clem: (future stages) stirring -> reaction on the diffused line
  // are composed here

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
        // Fernando-Clem: drift(e) now reflects the PHYSICAL p dV work of the
        // isobaric line (rho breathes), not a numerical error - drift(h) above
        // is the conservation check for the constant-pressure T-transport
        amrex::Print() << ", drift(e) = " << diag.diffusion.drift_energy
                       << ", drift(Y) = " << diag.diffusion.drift_species;
      }
      amrex::Print() << '\n';
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

// ================================================================================
// END OF CONSOLIDATED BUNDLE
// ================================================================================
