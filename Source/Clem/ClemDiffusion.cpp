#include <algorithm>
#include <cmath>
#include <limits>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Random.H>

#include "ClemDiffusion.H"
#include "ClemEosUtil.H"
#include "ClemStirring.H"

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

void
DiffusionOperator::tripletMap(
  LineScratch& lscr, const int start, const int n_seg)
{
  const int n = lscr.n_elem;
  if (n_seg < 6 || n_seg > n || (n_seg % 3) != 0) {
    return; // Fernando-Clem: not a mappable segment - eddyElements guards this
  }

  // Fernando-Clem: the permutation, as an offset table within the segment.
  // Three compressed copies of the segment, the middle one reversed:
  //   first  third -> offsets 0, 3, 6, ...
  //   middle third -> offsets n_seg-2, n_seg-5, ...   (reversed)
  //   last   third -> offsets 2, 5, 8, ...
  // Every offset in [0, n_seg) appears exactly once (verified exhaustively for
  // every legal n_seg and start), which is why the map moves fluid without
  // moving mass on an equal-mass line and contributes no conservation error of
  // its own - see the header
  const int n3 = n_seg / 3;
  static thread_local amrex::Vector<int> src;
  static thread_local amrex::Vector<amrex::Real> buf;
  src.resize(n_seg);
  buf.resize(static_cast<std::size_t>(n_seg) * (4 + NUM_SPECIES));

  for (int j = 0; j < n3; ++j) {
    src[j] = 3 * j;
    src[n3 + j] = n_seg - 2 - 3 * j;
    src[2 * n3 + j] = 3 * j + 2;
  }

  // Fernando-Clem: read the whole segment out first - the map is not in-place
  // (a destination can be read after it has been written otherwise)
  const int stride = 4 + NUM_SPECIES;
  for (int m = 0; m < n_seg; ++m) {
    const int l = (start + src[m]) % n; // periodic line: the segment wraps
    amrex::Real* b = &buf[static_cast<std::size_t>(m) * stride];
    b[0] = lscr.rho[l];
    b[1] = lscr.T[l];
    b[2] = lscr.eint[l];
    b[3] = lscr.press[l];
    for (int k = 0; k < NUM_SPECIES; ++k) {
      b[4 + k] = lscr.Y[l * NUM_SPECIES + k];
    }
  }

  // Fernando-Clem: element MASS is deliberately not permuted. The regrid leaves
  // the line equal-mass, so the elements the map exchanges hold identical dm
  // and permuting the fluid state alone is the mass-coordinate triplet map
  for (int m = 0; m < n_seg; ++m) {
    const int l = (start + m) % n;
    const amrex::Real* b = &buf[static_cast<std::size_t>(m) * stride];
    lscr.rho[l] = b[0];
    lscr.T[l] = b[1];
    lscr.eint[l] = b[2];
    lscr.press[l] = b[3];
    for (int k = 0; k < NUM_SPECIES; ++k) {
      lscr.Y[l * NUM_SPECIES + k] = b[4 + k];
    }
  }
}

DiffusionOperator::Diagnostics
DiffusionOperator::diffuse(
  ClemManager& mgr,
  const int lev,
  const amrex::Real dt,
  const TransParmType* tparm,
  ReactionOperator::Diagnostics* react_diag,
  const amrex::MultiFab* nut_strain,
  const amrex::Real filter_width)
{
  BL_PROFILE("clem::DiffusionOperator::diffuse()");
  AMREX_ALWAYS_ASSERT(mgr.isDefined());
  AMREX_ALWAYS_ASSERT(tparm != nullptr);

  // Fernando-Clem: stirring is on only if the caller built the (nu_t, |S|)
  // field for it. Delta must then be a real length
  const bool do_stir = (nut_strain != nullptr) && (filter_width > 0.0);
  AMREX_ALWAYS_ASSERT( !do_stir || nut_strain->nComp() >= NutComp::ncomp);

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
    // Fernando-Clem: the cell's LES turbulence, built once per step by
    // clem::computeTurbViscosity (same BoxArray as the state). Held BY VALUE -
    // an empty Array4 when stirring is off, and never dereferenced then
    const amrex::Array4<const amrex::Real> nut_arr =
      do_stir ? nut_strain->const_array(grid_id)
              : amrex::Array4<const amrex::Real>{};

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      const auto& plist = mgr.particlesInCell(lev, grid_id, iv);
      if (!gatherLine(lscr, plist, particles)) {
        ++diag.n_skipped;
        continue;
      }

      // Fernando-Clem: line enthalpy before diffusion, h = e + p/rho (raw
      // gathered state - does not depend on updateProperties)
      for (int l = 0; l < lscr.n_elem; ++l) {
        H_before += lscr.mass[l] * (lscr.eint[l] + lscr.press[l] / lscr.rho[l]);
      }

      // Fernando-Clem: ADAPTIVE sub-cycle - Maxwell [M §3.3.3] steps 14-20.
      // Delta t_diff is recomputed FRESH at the top of every iteration, from
      // whatever (T, Y) the line is currently at (post-diffusion AND, when
      // react_diag != nullptr, post-reaction too), then clipped to whatever
      // remains of dt. This is NOT "compute Delta t_diff once, divide dt into
      // N equal pieces": reaction can release heat fast enough mid-line to
      // change the diffusion stability limit substep to substep, and a
      // stale, pre-reaction estimate would silently under- or over-resolve
      // later substeps. The loop always covers dt exactly (or is force-
      // completed and flagged once n_capped substeps is reached - see below)
      amrex::Real t_clem      = 0.0;
      int n_sub               = 0;
      bool capped_this_line   = false;
      // Fernando-Clem: co-stepping bookkeeping - a line counts as "retried"
      // (react_diag->n_retried_lines) if ANY of its substeps needed a
      // reactClemLine retry, not once per substep
      bool line_retried               = false;
      constexpr amrex::Real time_tol  = 1.0e-12;
      const amrex::Real tol           = time_tol * dt;

      // Fernando-Clem: ===================== STIRRING =====================
      // The stirring clock, set up ONCE per line and per LES step. dt_stir is
      // LATCHED here and never recomputed inside the sub-cycle: its inputs are
      // the cell's nu_t and |S|, which are properties of the RESOLVED field and
      // are frozen while the subgrid sub-cycles. Only Delta t_diff changes from
      // substep to substep, because only the transport coefficients do.
      //
      //   t_stir = t + Delta t_stir   is the next epoch
      //   eta                          feeds the eddy-size PDF at each epoch
      //
      // An inactive closure (no subgrid turbulence, or the cell is resolved so
      // there is no inertial range) leaves dt_stir = +max, the first epoch
      // never arrives, and the loop below is EXACTLY the pre-stirring one
      amrex::Real dt_stir = std::numeric_limits<amrex::Real>::max();
      amrex::Real t_stir  = std::numeric_limits<amrex::Real>::max();
      amrex::Real eta     = 0.0;
      int n_stir          = 0;
      if (do_stir) {
        const amrex::Real nu   = lineViscosity(lscr, tparm);
        const amrex::Real nu_t = nut_arr(iv, NutComp::nut);
        const amrex::Real smag = nut_arr(iv, NutComp::strain);
        dt_stir = computeDtStirring(nu, filter_width, nu_t, smag);
        eta     = kolmogorovScale(nu, nu_t, smag);
        t_stir  = dt_stir; // first epoch: t = 0 at the start of the sub-cycle
      }

      while (dt - t_clem > time_tol * dt) {
        // Fernando-Clem: coefficients AND the stability limit are re-derived
        // from the current (T, Y) every iteration - both must be fresh
        // together, since stableSubStep reads updateProperties' output
        updateProperties(lscr, tparm);

        amrex::Real dt_sub              = stableSubStep(lscr, area2);
        const amrex::Real remaining     = dt - t_clem;
        bool forced_to_end              = false;
        if (dt_sub >= remaining) {
          // Fernando-Clem: diffusion is not (or no longer) the restrictive
          // constraint on what is left of dt - finish the line in this step
          dt_sub = remaining;
        } else if (n_sub + 1 >= max_sub) {
          // Fernando-Clem: about to need MORE than clem.diffusion_max_substeps
          // iterations to cover dt at the true CFL limit. Force completion in
          // this (oversized, above-CFL) step rather than looping forever -
          // silently capping is how an unstable line turns into a NaN three
          // steps later, so it is counted and shouted instead
          dt_sub            = remaining;
          capped_this_line  = true;
          forced_to_end     = true;
        }

        // Fernando-Clem: [stir] the THIRD clip - to the next stirring epoch.
        // dt_sub is now min(Delta t_diff, what remains of dt, t_stir - t).
        // Clipping here is what lands the map exactly ON its epoch instead of
        // smearing it into the middle of a transport step, and it is also what
        // guarantees at most ONE map per substep. Skipped when the substep cap
        // has already forced this step to the end of dt - that step is a
        // damage-limitation step, and shortening it again would not terminate
        if (!forced_to_end && t_clem + dt_sub > t_stir) {
          dt_sub = t_stir - t_clem;
        }

        computeRates(lscr, area2);
        advanceLine(lscr, dt_sub, diag);

        // Fernando-Clem: reaction co-stepped RIGHT AFTER diffusion, on the
        // SAME dt_sub (Maxwell [M §3.3.3] steps 14-20) - not deferred to a
        // separate stage over the full dt. Operates on lscr directly; the
        // particle write-back happens once, below, after the whole sub-cycle
        // (diffusion+reaction together) has covered dt
        if (react_diag != nullptr) {
          const int n_retries = ReactionOperator::reactSubstep(lscr, dt_sub, *react_diag);
          if (n_retries > 0) {
            line_retried = true;
          }
        }

        t_clem += dt_sub;
        ++n_sub;

        // Fernando-Clem: [stir] the epoch has been reached - fire ONE triplet
        // map, then push the epoch forward by the (constant) Delta t_stir.
        // Order matters: the map comes AFTER diffusion+reaction have advanced
        // the line to t_stir, so it rearranges the state the eddy would
        // actually have found there
        if (t_clem >= t_stir - tol) {
          const amrex::Real l = sampleEddySize(eta, filter_width, amrex::Random());
          const int n_seg     = eddyElements(l, filter_width, lscr.n_elem);
          if (n_seg > 0) {
            // Fernando-Clem: the sampler gives the CENTRE; the map wants the
            // first element. The line is periodic, so this wraps
            const int centre  = sampleEddyCenter(amrex::Random(), lscr.n_elem);
            const int start   = ((centre - n_seg / 2) % lscr.n_elem + lscr.n_elem) % lscr.n_elem;
            tripletMap(lscr, start, n_seg);
            ++n_stir;
          } else {
            // Fernando-Clem: the eddy occurred but the line is too coarse to
            // represent it (< 6 elements). It is REJECTED, not rounded up to 6
            // - rounding up would inject mixing at a scale the closure never
            // asked for. The clock still advances: the event happened
            ++diag.n_stir_rejected;
          }
          t_stir += dt_stir;
        }
      }

      if (capped_this_line) {
        ++diag.n_capped;
      }
      // Fernando-Clem: [stir] the substep cap jumped the sub-cycle to the end
      // of dt with epochs still pending - those maps were silently dropped, so
      // the realised D_T on this line is below target for reasons that have
      // nothing to do with the closure
      if (capped_this_line && do_stir && t_stir < dt) {
        ++diag.n_stir_truncated;
      }
      diag.total_stir_events += static_cast<amrex::Long>(n_stir);
      diag.max_stir_events = amrex::max<int>(diag.max_stir_events, n_stir);
      diag.max_substeps = amrex::max<int>(diag.max_substeps, n_sub);
      // Fernando-Clem: accumulate the TOTAL work (summed across ranks below):
      // one unit = one full-line computeRates+advanceLine sweep
      diag.total_substeps += static_cast<amrex::Long>(n_sub);
      ++diag.n_lines;

      if (react_diag != nullptr) {
        ++react_diag->n_lines;
        if (line_retried) {
          ++react_diag->n_retried_lines;
        }
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
  amrex::ParallelDescriptor::ReduceLongSum(diag.total_stir_events);
  amrex::ParallelDescriptor::ReduceIntMax(diag.max_stir_events);
  amrex::ParallelDescriptor::ReduceLongSum(diag.n_stir_rejected);
  amrex::ParallelDescriptor::ReduceLongSum(diag.n_stir_truncated);

  // Fernando-Clem: co-stepped reaction diagnostics - diffuse() is the only
  // place reactSubstep() is called, so it owns this reduction (react()'s
  // standalone path does the same reduction itself, at the end of react())
  if (react_diag != nullptr) {
    amrex::ParallelDescriptor::ReduceLongSum(react_diag->n_lines);
    amrex::ParallelDescriptor::ReduceRealMax(react_diag->max_neg_Y);
    amrex::ParallelDescriptor::ReduceIntMax(react_diag->max_retries);
    amrex::ParallelDescriptor::ReduceLongSum(react_diag->n_retried_lines);
  }

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
