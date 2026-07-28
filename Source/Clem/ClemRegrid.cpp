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
