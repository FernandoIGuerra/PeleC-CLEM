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
