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
