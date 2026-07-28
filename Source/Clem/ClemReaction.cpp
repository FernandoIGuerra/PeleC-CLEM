#include <algorithm>
#include <vector>

#include <AMReX_ParallelDescriptor.H>

#include "ReactorCvode.H"
#include "ReactorTypes.H"

#include "ClemReaction.H"

namespace clem {

namespace {
// Fernando-Clem: the persistent reactor's CVODE vector is sized once, at
// initClemReactor(config::n_lem) time, from PelePhysics's NUM_LEM (see
// ReactorClemConfig.H). This is the one place both constants are visible
// together - a silent mismatch would size the batched vector wrong and
// corrupt every line, so it is a hard compile error, not a runtime check.
static_assert(
  NUM_LEM == config::n_lem,
  "PelePhysics ReactorClemConfig.H's NUM_LEM must equal "
  "clem::config::n_lem (Source/Clem/ClemIndex.H) - rebuild with matching "
  "-DNUM_LEM=<n> and -DCLEM_NLEM=<n>.");
} // namespace

bool
ReactionOperator::gatherLine(
  const std::vector<int>& plist,
  ClemParticleContainer::ParticleType* particles,
  ReactionScratch& scr)
{
  const int n = static_cast<int>(plist.size());
  if (n < config::n_lem) {
    return false; // under-populated cell (shouldn't happen post-regrid)
  }

  // Fernando-Clem: reused scratch (like LineScratch), sorted into line order
  static thread_local std::vector<int> sorted;
  sorted.assign(plist.begin(), plist.end());
  std::sort(sorted.begin(), sorted.end(), [&](const int a, const int b) {
    return particles[a].idata(IntData::lem_index) <
           particles[b].idata(IntData::lem_index);
  });

  auto eos = pele::physics::PhysicsType::eos();

  for (int l = 0; l < config::n_lem; ++l) {
    const auto& p = particles[sorted[l]];
    const amrex::Real rho = p.rdata(RealData::rho);
    const amrex::Real T = p.rdata(RealData::T);
    const amrex::Real p_held = p.rdata(RealData::press);
    const amrex::Real m = p.rdata(RealData::mass);
    // Fernando-Clem: same placeholder guard as DiffusionOperator::gatherLine
    // - a padded/empty-cell element would poison the whole batched line
    if (m <= 0.0 || rho <= 0.0 || p_held <= 0.0) {
      return false;
    }

    amrex::Real Y[NUM_SPECIES];
    for (int sp = 0; sp < NUM_SPECIES; ++sp) {
      Y[sp] = p.rdata(RealData::Y0 + sp);
    }

    // Fernando-Clem: HP mode - conservatives[][NUM_SPECIES] must carry rho*h
    // (specific enthalpy density), not rho*e, so the reactor's RHY2T
    // inversion recovers T consistently with the element's HELD pressure
    amrex::Real e = 0.0;
    eos.RTY2E(rho, T, Y, e);
    const amrex::Real h = e + p_held / rho;

    scr.pidx[l] = sorted[l];
    for (int sp = 0; sp < NUM_SPECIES; ++sp) {
      scr.conservatives[l][sp] = rho * Y[sp];
      scr.rhs[l][sp] = 0.0; // sequential splitting - no external source
    }
    scr.conservatives[l][NUM_SPECIES] = rho * h;
    scr.rhs[l][NUM_SPECIES] = 0.0;
    scr.primitives[l][NUM_SPECIES] = T; // Newton warm-start only (see flatten)
    scr.density[l] = rho;
    scr.press_held[l] = p_held;
    scr.mass[l] = m;
  }
  return true;
}

void
ReactionOperator::scatterLine(
  ClemParticleContainer::ParticleType* particles,
  const ReactionScratch& scr,
  Diagnostics& diag)
{
  auto eos = pele::physics::PhysicsType::eos();

  for (int l = 0; l < config::n_lem; ++l) {
    auto& p = particles[scr.pidx[l]];

    // Fernando-Clem: post-reaction composition. Reaction conserves total
    // mass exactly, so rho_sum should already equal the pre-reaction rho to
    // solver tolerance; recomputed fresh here (not assumed) - same
    // discipline as DiffusionOperator's rho recovery
    amrex::Real rho_sum = 0.0;
    for (int sp = 0; sp < NUM_SPECIES; ++sp) {
      rho_sum += scr.conservatives[l][sp];
    }
    const amrex::Real one_rho = (rho_sum > 0.0) ? 1.0 / rho_sum : 0.0;

    amrex::Real Y[NUM_SPECIES];
    amrex::Real worst_neg = 0.0;
    for (int sp = 0; sp < NUM_SPECIES; ++sp) {
      amrex::Real Yk = scr.conservatives[l][sp] * one_rho;
      if (Yk < 0.0) {
        worst_neg = amrex::max<amrex::Real>(worst_neg, -Yk);
        Yk = 0.0;
      }
      Y[sp] = Yk;
    }
    diag.max_neg_Y = amrex::max<amrex::Real>(diag.max_neg_Y, worst_neg);

    const amrex::Real T_new = scr.primitives[l][NUM_SPECIES];

    // Fernando-Clem: ISOBARIC closure - (rho, e) reconstructed at the HELD
    // element pressure from the reacted (T, Y), exactly like
    // DiffusionOperator::advanceLine. Pressure itself is left unchanged (it
    // is still the reference p_1 the next isentropic LES->SGS stage
    // compresses from).
    amrex::Real rho_new = 0.0;
    amrex::Real e_new = 0.0;
    eos.PYT2RE(scr.press_held[l], Y, T_new, rho_new, e_new);

    p.rdata(RealData::T) = T_new;
    p.rdata(RealData::rho) = rho_new;
    p.rdata(RealData::eint) = e_new;
    p.rdata(RealData::vol) = scr.mass[l] / rho_new;
    for (int sp = 0; sp < NUM_SPECIES; ++sp) {
      p.rdata(RealData::Y0 + sp) = Y[sp];
    }
  }
}

pele::physics::reactions::ReactorCvode&
ReactionOperator::reactor()
{
  // Fernando-Clem: ONE persistent reactor per rank (matches LineScratch's
  // per-rank scratch-arena pattern), shared by BOTH react() and
  // reactSubstep(). init()/initClemReactor() run ONCE, lazily, on whichever
  // entry point is called first; reused via CVodeReInit for every subsequent
  // line for the rest of the run. NOT thread_local: CLEM's particle loops
  // are not currently run inside an OpenMP region (plain MPI-rank
  // iteration) - if that changes, this must become thread_local like
  // LineScratch/ReactionScratch already are.
  static pele::physics::reactions::ReactorCvode instance;
  static bool initialized = false;
  if (!initialized) {
    instance.init(  pele::physics::reactions::ReactorTypes::h_reactor_type, config::n_lem);
    instance.initClemReactor(config::n_lem);
    initialized = true;
  }
  return instance;
}

int
ReactionOperator::reactSubstep(
  LineScratch& lscr, const amrex::Real dt_sub, Diagnostics& diag)
{
  AMREX_ASSERT(lscr.n_elem == config::n_lem);

  auto eos = pele::physics::PhysicsType::eos();
  auto& scr = ReactionScratch::get();

  // Fernando-Clem: gather FROM the LineScratch diffusion already populated -
  // same primitive->conservative(HP) packing as gatherLine, just sourced from
  // lscr instead of the particle AoS
  for (int l = 0; l < config::n_lem; ++l) {
    const amrex::Real rho = lscr.rho[l];
    const amrex::Real T = lscr.T[l];
    const amrex::Real p_held = lscr.press[l];
    const amrex::Real m = lscr.mass[l];

    amrex::Real Y[NUM_SPECIES];
    for (int sp = 0; sp < NUM_SPECIES; ++sp) {
      Y[sp] = lscr.Y[l * NUM_SPECIES + sp];
    }

    amrex::Real e = 0.0;
    eos.RTY2E(rho, T, Y, e);
    const amrex::Real h = e + p_held / rho;

    for (int sp = 0; sp < NUM_SPECIES; ++sp) {
      scr.conservatives[l][sp] = rho * Y[sp];
      scr.rhs[l][sp] = 0.0; // sequential splitting - no external source
    }
    scr.conservatives[l][NUM_SPECIES] = rho * h;
    scr.rhs[l][NUM_SPECIES] = 0.0;
    scr.primitives[l][NUM_SPECIES] = T; // Newton warm-start only
    scr.density[l] = rho;
    scr.press_held[l] = p_held;
    scr.mass[l] = m;
  }

  amrex::Real dt_react = dt_sub;
  amrex::Real time = 0.0;
  const int n_retries = reactor().reactClemLine(
    scr.primitives, scr.conservatives, scr.rhs, scr.density, dt_react, time);
  diag.max_retries = amrex::max(diag.max_retries, n_retries);

  // Fernando-Clem: scatter BACK into the LineScratch (not the particles) -
  // diffusion's own scatterLine, called once after ALL substeps for this
  // line, writes the final (T, rho, eint, Y, vol) to the particles
  for (int l = 0; l < config::n_lem; ++l) {
    amrex::Real rho_sum = 0.0;
    for (int sp = 0; sp < NUM_SPECIES; ++sp) {
      rho_sum += scr.conservatives[l][sp];
    }
    const amrex::Real one_rho = (rho_sum > 0.0) ? 1.0 / rho_sum : 0.0;

    amrex::Real Y[NUM_SPECIES];
    amrex::Real worst_neg = 0.0;
    for (int sp = 0; sp < NUM_SPECIES; ++sp) {
      amrex::Real Yk = scr.conservatives[l][sp] * one_rho;
      if (Yk < 0.0) {
        worst_neg = amrex::max<amrex::Real>(worst_neg, -Yk);
        Yk = 0.0;
      }
      Y[sp] = Yk;
    }
    diag.max_neg_Y = amrex::max<amrex::Real>(diag.max_neg_Y, worst_neg);

    const amrex::Real T_new = scr.primitives[l][NUM_SPECIES];
    amrex::Real rho_new = 0.0;
    amrex::Real e_new = 0.0;
    eos.PYT2RE(scr.press_held[l], Y, T_new, rho_new, e_new);

    lscr.T[l] = T_new;
    lscr.rho[l] = rho_new;
    lscr.eint[l] = e_new;
    for (int sp = 0; sp < NUM_SPECIES; ++sp) {
      lscr.Y[l * NUM_SPECIES + sp] = Y[sp];
    }
  }
  return n_retries;
}

ReactionOperator::Diagnostics
ReactionOperator::react(ClemManager& mgr, const int lev, const amrex::Real dt)
{
  BL_PROFILE("clem::ReactionOperator::react()");
  AMREX_ALWAYS_ASSERT(mgr.isDefined());

  auto& pc = mgr.particleContainer();
  if (!mgr.IsCellMapValid()) {
    mgr.rebuildCellMap(lev);
  }

  Diagnostics diag;
  // Fernando-Clem: per-rank persistent arena (ReactionScratch::get(), see
  // ClemReaction.H) - no per-call construction/zero-init cost, same
  // discipline as ClemDiffusion's LineScratch
  auto& scr = ReactionScratch::get();

  for (ClemParIter pti(pc, lev); pti.isValid(); ++pti) {
    const int grid_id         = pti.index();
    const amrex::Box& bx      = pti.validbox();
    auto* particles           = pti.GetArrayOfStructs().data();

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      const auto& plist = mgr.particlesInCell(lev, grid_id, iv);
      if (!gatherLine(plist, particles, scr)) {
        ++diag.n_skipped;
        continue;
      }

      amrex::Real dt_react    = dt;
      amrex::Real time        = 0.0;
      const int n_retries     = reactor().reactClemLine(scr.primitives, scr.conservatives, scr.rhs, scr.density, dt_react,time);
      
      if (n_retries > 0) {
        diag.max_retries = amrex::max(diag.max_retries, n_retries);
        ++diag.n_retried_lines;
      }

      scatterLine(particles, scr, diag);
      ++diag.n_lines;
    }
  }

  amrex::ParallelDescriptor::ReduceLongSum(diag.n_lines);
  amrex::ParallelDescriptor::ReduceLongSum(diag.n_skipped);
  amrex::ParallelDescriptor::ReduceRealMax(diag.max_neg_Y);
  amrex::ParallelDescriptor::ReduceIntMax(diag.max_retries);
  amrex::ParallelDescriptor::ReduceLongSum(diag.n_retried_lines);
  return diag;
}

} // namespace clem
