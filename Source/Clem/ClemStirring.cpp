#include <cmath>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include "ClemManager.H"
#include "ClemStirring.H"

namespace clem {

StirParams StirParams::m_params;

void
StirParams::readParams()
{
  amrex::ParmParse pp("clem");
  pp.query("do_stir", m_params.do_stir);

  // Fernando-Clem: a triplet map is instantaneous, so it does not get a stage
  // of its own - it rides inside the diffusion sub-cycle as a clip on
  // Delta t_diff (DiffusionOperator::diffuse). With clem.do_diffusion = 0 there
  // is no sub-cycle to carry the clock and stirring would be silently inert,
  // which is worse than refusing to start
  if (m_params.do_stir && !ClemManager::doDiffusion()) {
    amrex::Abort(
      "clem.do_stir = 1 requires clem.do_diffusion = 1: the stirring clock is "
      "carried by the diffusion sub-cycle (Delta t_diff is clipped to the next "
      "stirring epoch), so with diffusion off no triplet map would ever fire.");
  }

  if (m_params.do_stir && ClemManager::verbose() > 0) {
    amrex::Print()
      << "CLEM stirring: ON - D_T = nu_t and eta = (nu^3/(nu_t |S|^2))^(1/4), "
         "both from PeleC's active LES model (Sc_t = 1)\n";
  }
}

void
checkLesCompatibility(const LesView& les)
{
  if (!StirParams::get().do_stir) {
    return;
  }
  // Fernando-Clem: what stirring actually needs is an eddy-viscosity
  // COEFFICIENT, not PeleC's LES term. computeTurbViscosity forms nu_t itself,
  // from that coefficient and the resolved velocity gradients, and never calls
  // into PeleC::getLESTerm - so the requirement here is on C_s (or C_w), not on
  // pelec.do_les.
  //
  // That distinction is not pedantry, it is what makes the stage usable:
  // PeleC's own LES stage aborts for NUM_SPECIES > 2 (PeleC.cpp, init_les),
  // so demanding pelec.do_les = 1 would make clem.do_stir impossible for every
  // real mechanism. The modelling argument for wanting a subgrid MOMENTUM
  // closure alongside this subgrid SCALAR closure still stands - if PeleC's LES
  // ever supports multi-component systems, turn pelec.do_les on as well - but
  // it cannot be enforced here without disabling the feature outright.
  if (les.model < 0 || les.model > 3) {
    amrex::Abort("clem stirring: unrecognised LES model");
  }
  if (les.model == 3) {
    amrex::Abort(
      "clem.do_stir does not implement a cell-centred nu_t for les_model = 3 "
      "(Vreman). Use les_model 0 (Smagorinsky), 1 (dynamic Smagorinsky) or 2 "
      "(WALE).");
  }
  if (les.model == 1 && les.coeffs == nullptr) {
    // Fernando-Clem: the dynamic model is the ONE case that genuinely needs
    // PeleC's LES machinery - the per-cell C_s^2 field is built by init_les()
    amrex::Abort(
      "clem stirring: les_model = 1 (dynamic Smagorinsky) needs PeleC's "
      "per-cell C_s^2 field, which only exists when pelec.do_les = 1. Set "
      "pelec.do_les = 1, or use les_model 0 (constant C_s = pelec.Cs), which "
      "CLEM evaluates on its own.");
  }
  // Fernando-Clem: a zero coefficient means nu_t == 0 everywhere, i.e. the
  // closure is inert and no map ever fires. Caught here rather than discovered
  // as "stirring is on and nothing happens" after a long run
  if (les.model == 0 && les.Cs <= 0.0) {
    amrex::Abort(
      "clem.do_stir with les_model = 0 needs pelec.Cs > 0: nu_t = C_s^2 "
      "Delta^2 |S| would be identically zero and no triplet map would fire.");
  }
  if (les.model == 2 && les.Cw <= 0.0) {
    amrex::Abort(
      "clem.do_stir with les_model = 2 (WALE) needs pelec.Cw > 0: nu_t would "
      "be identically zero and no triplet map would fire.");
  }
}

amrex::Real
lineViscosity(const LineScratch& lscr, const TransParmType* tparm)
{
  AMREX_ALWAYS_ASSERT(tparm != nullptr);
  auto trans = pele::physics::PhysicsType::transport();

  // Fernando-Clem: mu only - lambda and Ddiag are the diffusion operator's
  // business and are already recomputed every sub-step there
  constexpr bool get_xi = false;
  constexpr bool get_mu = true;
  constexpr bool get_lam = false;
  constexpr bool get_Ddiag = false;
  constexpr bool get_chi = false;

  amrex::Real m_tot = 0.0;
  amrex::Real m_nu = 0.0;
  for (int l = 0; l < lscr.n_elem; ++l) {
    const amrex::Real rho_l = lscr.rho[l];
    if (rho_l <= 0.0) {
      continue;
    }
    amrex::Real Yl[NUM_SPECIES];
    for (int k = 0; k < NUM_SPECIES; ++k) {
      Yl[k] = lscr.Y[l * NUM_SPECIES + k];
    }

    amrex::Real mu = 0.0;
    amrex::Real xi = 0.0;
    amrex::Real lam = 0.0;
    amrex::Real Dl[NUM_SPECIES] = {0.0};
    amrex::Real* chi_mix = nullptr;
    trans.transport(
      get_xi, get_mu, get_lam, get_Ddiag, get_chi, lscr.T[l], rho_l, Yl, Dl,
      chi_mix, mu, xi, lam, tparm);

    // Fernando-Clem: MASS-weighted, because the LEM line is a mass grid - a
    // volume average would over-weight the hot, light elements
    m_nu += lscr.mass[l] * mu / rho_l;
    m_tot += lscr.mass[l];
  }
  return (m_tot > 0.0) ? m_nu / m_tot : 0.0;
}

amrex::Real
filterWidth(const amrex::Geometry& geom)
{
  // Fernando-Clem: Delta = (dx dy dz)^(1/3), geometric mean over the
  // dimensions that exist
  const auto* dx = geom.CellSize();
  amrex::Real vol = 1.0;
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    vol *= dx[d];
  }
  return std::pow(vol, 1.0 / static_cast<amrex::Real>(AMREX_SPACEDIM));
}

void
computeTurbViscosity(
  const amrex::MultiFab& state,
  const int rho_indx,
  const int mom_indx,
  const amrex::Geometry& geom,
  const LesView& les,
  amrex::MultiFab& nut_strain)
{
  BL_PROFILE("clem::computeTurbViscosity()");
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    state.nGrow() >= 1,
    "clem::computeTurbViscosity needs >= 1 ghost cell on the state (central "
    "differences of the resolved velocity)");
  AMREX_ALWAYS_ASSERT(nut_strain.nComp() >= NutComp::ncomp);

  const auto* dx = geom.CellSize();
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> half_odx{};
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    half_odx[d] = 0.5 / dx[d]; // central difference: (f(+1)-f(-1))/(2 dx)
  }

  const amrex::Real delta = filterWidth(geom);
  const amrex::Real delta2 = delta * delta;
  const int model = les.model;
  // Fernando-Clem: models 0 and 1 differ ONLY in where C_s^2 comes from, so the
  // constant is hoisted here and the per-cell field read inside the loop
  const amrex::Real cs2_const = les.Cs * les.Cs;
  const amrex::Real cw2 = les.Cw * les.Cw;
  const bool dynamic = (model == 1);
  const int cs2_comp = les.cs2_comp;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(nut_strain, amrex::TilingIfNotGPU()); mfi.isValid();
       ++mfi) {
    const amrex::Box& bx = mfi.tilebox();
    const auto& s = state.const_array(mfi);
    const auto& out = nut_strain.array(mfi);
    // Fernando-Clem: dummy array when the model is not dynamic - Array4 cannot
    // be conditionally constructed inside the device lambda
    const auto& cs2f = dynamic ? les.coeffs->const_array(mfi) : s;

    amrex::ParallelFor(
      bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};

        // Fernando-Clem: g[m][n] = d u_m / d x_n. Only the components that
        // exist in this dimensionality are filled; the rest stay zero, which is
        // what the 2-D/1-D invariants need
        amrex::Real g[AMREX_SPACEDIM][AMREX_SPACEDIM] = {{0.0}};
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
          const auto ivp = iv + amrex::IntVect::TheDimensionVector(n);
          const auto ivm = iv - amrex::IntVect::TheDimensionVector(n);
          const amrex::Real rp = s(ivp, rho_indx);
          const amrex::Real rm = s(ivm, rho_indx);
          if (rp <= 0.0 || rm <= 0.0) {
            continue; // Fernando-Clem: unusable neighbour, leave column zero
          }
          for (int m = 0; m < AMREX_SPACEDIM; ++m) {
            const amrex::Real up = s(ivp, mom_indx + m) / rp;
            const amrex::Real um = s(ivm, mom_indx + m) / rm;
            g[m][n] = (up - um) * half_odx[n];
          }
        }

        // Fernando-Clem: S_ij = (1/2)(g_ij + g_ji), |S| = sqrt(2 S_ij S_ij)
        amrex::Real sij_sij = 0.0;
        for (int m = 0; m < AMREX_SPACEDIM; ++m) {
          for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            const amrex::Real sij = 0.5 * (g[m][n] + g[n][m]);
            sij_sij += sij * sij;
          }
        }
        const amrex::Real smag = std::sqrt(2.0 * sij_sij);

        amrex::Real nut = 0.0;
        if (model == 2) {
          // Fernando-Clem: WALE. g2_ij = g_ik g_kj; the traceless symmetric
          // part Sd_ij = (1/2)(g2_ij + g2_ji) - (1/3) delta_ij g2_kk
          amrex::Real g2[AMREX_SPACEDIM][AMREX_SPACEDIM] = {{0.0}};
          amrex::Real trace = 0.0;
          for (int m = 0; m < AMREX_SPACEDIM; ++m) {
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
              amrex::Real acc = 0.0;
              for (int q = 0; q < AMREX_SPACEDIM; ++q) {
                acc += g[m][q] * g[q][n];
              }
              g2[m][n] = acc;
            }
            trace += g2[m][m];
          }
          amrex::Real sd_sd = 0.0;
          for (int m = 0; m < AMREX_SPACEDIM; ++m) {
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
              const amrex::Real sd =
                0.5 * (g2[m][n] + g2[n][m]) - ((m == n) ? trace / 3.0 : 0.0);
              sd_sd += sd * sd;
            }
          }
          const amrex::Real den = std::pow(sij_sij, 2.5) + std::pow(sd_sd, 1.25);
          nut = (den > 0.0) ? cw2 * delta2 * std::pow(sd_sd, 1.5) / den : 0.0;
        } else {
          // Fernando-Clem: Smagorinsky (0) and dynamic Smagorinsky (1) -
          // nu_t = C_s^2 Delta^2 |S|, the ONLY difference being whether C_s^2
          // is the input constant or the per-cell filtered field
          const amrex::Real cs2 =
            dynamic ? amrex::max<amrex::Real>(cs2f(iv, cs2_comp), 0.0)
                    : cs2_const;
          nut = cs2 * delta2 * smag;
        }

        out(iv, NutComp::nut) = amrex::max<amrex::Real>(nut, 0.0);
        out(iv, NutComp::strain) = smag;
      });
  }
}

} // namespace clem
