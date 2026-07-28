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
bool ClemManager::m_do_react = false;

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
  pp.query("do_react", m_do_react);

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
                   << ", do_diffusion = " << m_do_diffusion
                   << ", do_react = " << m_do_react << '\n';
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
