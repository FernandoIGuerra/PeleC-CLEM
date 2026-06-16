
#include "ClemParticles/clem_diffusion.H"
#include <TransportParams.H>
#include "ClemParticles/clem_index_definition.H"
#include <limits>
#include <cmath>

// ============================================================
// clem_diffusion.cpp
//
// Each thermal diffusion scheme is an explicit full specialisation of
// diffusion::TemperatureDiffusionScheme<Scheme>.
//
// ALL schemes use CFL-based sub-stepping so they are unconditionally stable
// regardless of the LES timestep dt.  The difference between them is:
//   O2_Euler   — 2nd-order spatial stencil, Euler sub-steps  (r_max = 0.5)
//   O4_Euler   — 4th-order spatial stencil, Euler sub-steps  (r_max = 3/11)
//   O2_RK2     — 2nd-order spatial stencil, Heun  sub-steps  (r_max = 0.5)
//   O2_Substep — alias of O2_Euler (kept for backward compatibility)
//
// To add a new scheme:
//   1. Add an entry to DiffusionScheme in clem_numerics_config.H
//   2. Write a new template<> specialisation here
//   3. Select it in clem_setup.H
// ============================================================

namespace clem {

using setup::diffusion_scheme;

// -----------------------------------------------------------------------
// Per-thread scratch bank — allocated once per thread, reused every call.
// Avoids repeated stack allocation of fixed-size arrays across particles.
// thread_local is safe under both OpenMP (one thread per team member) and MPI.
// -----------------------------------------------------------------------
struct DiffScratch {
    diffusion::LemMatrix prim_work;
    diffusion::LemMatrix prim_star;     // O2_RK2 predictor stage
    diffusion::LemMatrix src_Y;         // species source rates  [NUM_LEM][sp]
    diffusion::LemMatrix src_Y_sum;     // accumulated species source
    diffusion::LemReal   src_T;         // thermal source at current sub-step
    diffusion::LemReal   src_T_sum;     // accumulated thermal source
    diffusion::LemReal   fluxR;         // total X-flux at right face (correction sum)
    diffusion::LemReal   fluxL;         // total X-flux at left face  (correction sum)
    diffusion::LemReal   src1;          // O2_RK2 stage 1
    diffusion::LemReal   src2;          // O2_RK2 stage 2
};
static thread_local DiffScratch tl_scratch;

// -----------------------------------------------------------------------
// Shared spatial-stencil helpers
// -----------------------------------------------------------------------

static void ComputeEnergySource_O2(
    const diffusion::LemReal&   density,
    const diffusion::LemMatrix& T_field,
    const diffusion::LemMatrix& transport,
    amrex::Real dm,
    diffusion::LemReal& EintSource)
{
    constexpr int N = NUM_LEM;
    for (int i = 0; i < N; i++) {
        const int iR = (i + 1) % N;
        const int iL = (i - 1 + N) % N;

        const amrex::Real rho  = density[i];
        const amrex::Real rhoR = density[iR];
        const amrex::Real rhoL = density[iL];

        const amrex::Real T  = T_field[i][NUM_SPECIES];
        const amrex::Real TR = T_field[iR][NUM_SPECIES];
        const amrex::Real TL = T_field[iL][NUM_SPECIES];

        const amrex::Real alphaR = diffusion::SafeHarmonicMean( rhoR * transport[iR][CLEM_dComp_lambda], rho * transport[i][CLEM_dComp_lambda]);
        const amrex::Real alphaL = diffusion::SafeHarmonicMean( rhoL * transport[iL][CLEM_dComp_lambda], rho * transport[i][CLEM_dComp_lambda]);
        EintSource[i] = (alphaR * (TR - T) - alphaL * (T - TL)) / dm;
    }
}

static void ComputeEnergySource_O4(
    const diffusion::LemReal&   density,
    const diffusion::LemMatrix& T_field,
    const diffusion::LemMatrix& transport,
    amrex::Real dm,
    diffusion::LemReal& EintSource)
{
    constexpr int N = NUM_LEM;
    constexpr amrex::Real c24 = 1.0 / 24.0;

    for (int i = 0; i < N; i++) {
        const int iR  = (i + 1) % N;
        const int iRR = (i + 2) % N;
        const int iL  = (i - 1 + N) % N;
        const int iLL = (i - 2 + N) % N;

        const amrex::Real rho = density[i];

        const amrex::Real T   = T_field[i][NUM_SPECIES];
        const amrex::Real TR  = T_field[iR][NUM_SPECIES];
        const amrex::Real TRR = T_field[iRR][NUM_SPECIES];
        const amrex::Real TL  = T_field[iL][NUM_SPECIES];
        const amrex::Real TLL = T_field[iLL][NUM_SPECIES];

        const amrex::Real alphaR = diffusion::SafeHarmonicMean( density[iR] * transport[iR][CLEM_dComp_lambda], rho * transport[i][CLEM_dComp_lambda]);
        const amrex::Real alphaL = diffusion::SafeHarmonicMean( density[iL] * transport[iL][CLEM_dComp_lambda], rho * transport[i][CLEM_dComp_lambda]);

        const amrex::Real F_R = alphaR * (-TRR + 27.0*TR - 27.0*T  + TL ) * c24 / dm;
        const amrex::Real F_L = alphaL * (-TR  + 27.0*T  - 27.0*TL + TLL) * c24 / dm;

        EintSource[i] = F_R - F_L;
    }
}

// Computes species diffusion sources for a given prim_work state.
// Outputs:
//   src_Y[i][sp]  — rate for Y[sp] at particle i  (written to rhs[i][sp])
//   src_T_sp[i]   — T rate contribution from enthalpy flux of species transport
// X (mole fractions) and H (enthalpies) are frozen at T^n (passed in as XH arrays).
static void ComputeSpeciesSource(
    const diffusion::LemReal&   density,
    const diffusion::LemMatrix& prim_work,
    const diffusion::LemMatrix& transport,
    amrex::Real dm,
    // frozen mole fractions and enthalpies computed once before sub-step loop
    const std::array<std::array<amrex::Real, NUM_SPECIES>, NUM_LEM>& X,
    const std::array<std::array<amrex::Real, NUM_SPECIES>, NUM_LEM>& H,
    diffusion::LemMatrix& src_Y,
    diffusion::LemReal&   src_T_sp,
    diffusion::LemReal&   fluxR,
    diffusion::LemReal&   fluxL)
{
    constexpr int N = NUM_LEM;

    // Reset outputs
    for (int i = 0; i < N; i++) {
        src_T_sp[i] = 0.0;
        fluxR[i]    = 0.0;
        fluxL[i]    = 0.0;
        for (int sp = 0; sp < NUM_SPECIES; sp++) src_Y[i][sp] = 0.0;
    }

    for (int i = 0; i < N; i++) {
        const int iR = (i + 1) % N;
        const int iL = (i - 1 + N) % N;
        const amrex::Real rho  = density[i];
        const amrex::Real rhoR = density[iR];
        const amrex::Real rhoL = density[iL];

        for (int sp = 0; sp < NUM_SPECIES; sp++) {
            const amrex::Real Xc = X[i][sp],  XR = X[iR][sp], XL = X[iL][sp];
            const amrex::Real Yc = prim_work[i][sp], YR = prim_work[iR][sp], YL = prim_work[iL][sp];

            const amrex::Real aR = diffusion::SafeHarmonicMean(rhoR*transport[iR][sp], rho*transport[i][sp]);
            const amrex::Real aL = diffusion::SafeHarmonicMean(rhoL*transport[iL][sp], rho*transport[i][sp]);

            const amrex::Real dX_R = (XR - Xc) / dm;
            const amrex::Real dX_L = (Xc - XL) / dm;

            fluxR[i] += aR * dX_R;
            fluxL[i] += aL * dX_L;

            src_Y[i][sp] = aR * dX_R - aL * dX_L;
            src_T_sp[i] += aR * dX_R * 0.25*(H[iR][sp]+H[i][sp])*(Yc+YR)
                         - aL * dX_L * 0.25*(H[iL][sp]+H[i][sp])*(Yc+YL);
        }
    }

    // Correction flux to enforce sum(mass flux) = 0
    for (int i = 0; i < N; i++) {
        const int iR = (i + 1) % N;
        const int iL = (i - 1 + N) % N;
        for (int sp = 0; sp < NUM_SPECIES; sp++) {
            const amrex::Real corr = 0.5*(prim_work[iR][sp]+prim_work[i][sp])*fluxR[i]
                                   - 0.5*(prim_work[iL][sp]+prim_work[i][sp])*fluxL[i];
            src_Y[i][sp] = src_Y[i][sp] - corr;
        }
    }
}

// r_max: scheme-dependent CFL coefficient.
// Stability condition for mass-coordinate stencil: dt <= r_max * dm² / (ρ · transport · A²)
// O2 schemes: r_max = 0.5; O4 scheme: r_max = 3/11.
static amrex::Real ComputeStableTimestep(
    const diffusion::LemReal&   density,
    const diffusion::LemMatrix& transport,
    amrex::Real dm,
    amrex::Real areaFactor,
    amrex::Real r_max = 0.5)
{
    const amrex::Real A2 = areaFactor * areaFactor;
    amrex::Real beta = 0.0;
    for (int i = 0; i < NUM_LEM; i++) {
        amrex::Real alpha_max = transport[i][CLEM_dComp_lambda];
        for (int sp = 0; sp < NUM_SPECIES; sp++)
            alpha_max = std::max(alpha_max, transport[i][sp]);
        beta = std::max(beta, density[i] * alpha_max * A2);
    }
    return (beta > 0.0) ? (r_max * dm * dm) / beta : std::numeric_limits<amrex::Real>::max();
}

// -----------------------------------------------------------------------
// Explicit template specialisations
// -----------------------------------------------------------------------

// ---- O2_Euler : 2nd-order space, Euler sub-steps (r_max = 0.5) --------
template<>
void diffusion::TemperatureDiffusionScheme<DiffusionScheme::O2_Euler>(
    LemReal& mass, LemReal& density,
    LemMatrix& primitives, LemMatrix& transport, LemMatrix& rhs,
    amrex::Real areaFactor, amrex::Real dt)
{
    constexpr int N = NUM_LEM;
    amrex::Real dm = 0.0;
    for (int i = 0; i < N; i++) dm += mass[i];
    dm /= static_cast<amrex::Real>(N);

    const amrex::Real scale  = areaFactor * areaFactor / dm;
    const amrex::Real dt_cfl = ComputeStableTimestep(density, transport, dm, areaFactor, 0.5);
    const int         N_sub  = std::max(1, static_cast<int>(std::ceil(dt / dt_cfl)));
    const amrex::Real dt_sub = dt / static_cast<amrex::Real>(N_sub);

    LemMatrix& prim_work = tl_scratch.prim_work;   prim_work = primitives;
    LemReal&   src_sum   = tl_scratch.src_T_sum;     src_sum   = {};
    LemReal&   src       = tl_scratch.src1;

    for (int sub = 0; sub < N_sub; sub++) {
        ComputeEnergySource_O2(density, prim_work, transport, dm, src);
        for (int i = 0; i < N; i++) {
            prim_work[i][NUM_SPECIES] += dt_sub * src[i] * scale;
            src_sum[i] += src[i];
        }
    }

    const amrex::Real inv_N = 1.0 / static_cast<amrex::Real>(N_sub);
    for (int i = 0; i < N; i++)
        rhs[i][NUM_SPECIES] += src_sum[i] * inv_N * scale;
}

// ---- O4_Euler : 4th-order space, Euler sub-steps (r_max = 3/11) -------
template<>
void diffusion::TemperatureDiffusionScheme<DiffusionScheme::O4_Euler>(
    LemReal& mass, LemReal& density,
    LemMatrix& primitives, LemMatrix& transport, LemMatrix& rhs,
    amrex::Real areaFactor, amrex::Real dt)
{
    constexpr int N = NUM_LEM;
    amrex::Real dm = 0.0;
    for (int i = 0; i < N; i++) dm += mass[i];
    dm /= static_cast<amrex::Real>(N);

    const amrex::Real scale  = areaFactor * areaFactor / dm;
    const amrex::Real r_max  = 3.0 / 11.0;
    const amrex::Real dt_cfl = ComputeStableTimestep(density, transport, dm, areaFactor, r_max);
    const int         N_sub  = std::max(1, static_cast<int>(std::ceil(dt / dt_cfl)));
    const amrex::Real dt_sub = dt / static_cast<amrex::Real>(N_sub);

    LemMatrix& prim_work = tl_scratch.prim_work;   prim_work = primitives;
    LemReal&   src_sum   = tl_scratch.src_T_sum;     src_sum   = {};
    LemReal&   src       = tl_scratch.src1;

    for (int sub = 0; sub < N_sub; sub++) {
        ComputeEnergySource_O4(density, prim_work, transport, dm, src);
        for (int i = 0; i < N; i++) {
            prim_work[i][NUM_SPECIES] += dt_sub * src[i] * scale;
            src_sum[i] += src[i];
        }
    }

    const amrex::Real inv_N = 1.0 / static_cast<amrex::Real>(N_sub);
    for (int i = 0; i < N; i++)
        rhs[i][NUM_SPECIES] += src_sum[i] * inv_N * scale;
}

// ---- O2_RK2 : 2nd-order space, Heun sub-steps (r_max = 0.5) -----------
// Each sub-step uses the predictor-corrector (Heun) method.
// Stability limit for RK2 on parabolic problems equals Euler's (r <= 0.5);
// the gain is 2nd-order temporal accuracy per sub-step.
template<>
void diffusion::TemperatureDiffusionScheme<DiffusionScheme::O2_RK2>(
    LemReal& mass, LemReal& density,
    LemMatrix& primitives, LemMatrix& transport, LemMatrix& rhs,
    amrex::Real areaFactor, amrex::Real dt)
{
    constexpr int N = NUM_LEM;
    amrex::Real dm = 0.0;
    for (int i = 0; i < N; i++) dm += mass[i];
    dm /= static_cast<amrex::Real>(N);

    const amrex::Real scale  = areaFactor * areaFactor / dm;
    const amrex::Real dt_cfl = ComputeStableTimestep(density, transport, dm, areaFactor, 0.5);
    const int         N_sub  = std::max(1, static_cast<int>(std::ceil(dt / dt_cfl)));
    const amrex::Real dt_sub = dt / static_cast<amrex::Real>(N_sub);

    LemMatrix& prim_work = tl_scratch.prim_work;   prim_work = primitives;
    LemMatrix& prim_star = tl_scratch.prim_star;
    LemReal&   src_sum   = tl_scratch.src_T_sum;     src_sum   = {};
    LemReal&   src1      = tl_scratch.src1;
    LemReal&   src2      = tl_scratch.src2;

    for (int sub = 0; sub < N_sub; sub++) {
        // Stage 1: stencil at T^(s)
        ComputeEnergySource_O2(density, prim_work, transport, dm, src1);

        // Predictor: T^* = T^(s) + dt_sub * src1 * scale
        prim_star = prim_work;
        for (int i = 0; i < N; i++)
            prim_star[i][NUM_SPECIES] += dt_sub * src1[i] * scale;

        // Stage 2: stencil at T^*
        ComputeEnergySource_O2(density, prim_star, transport, dm, src2);

        // Corrector: advance prim_work with averaged rate
        for (int i = 0; i < N; i++) {
            const amrex::Real src_avg = 0.5 * (src1[i] + src2[i]);
            prim_work[i][NUM_SPECIES] += dt_sub * src_avg * scale;
            src_sum[i] += src_avg;
        }
    }

    const amrex::Real inv_N = 1.0 / static_cast<amrex::Real>(N_sub);
    for (int i = 0; i < N; i++)
        rhs[i][NUM_SPECIES] += src_sum[i] * inv_N * scale;
}

// ---- O2_Substep : alias of O2_Euler (kept for backward compatibility) --
template<>
void diffusion::TemperatureDiffusionScheme<DiffusionScheme::O2_Substep>(
    LemReal& mass, LemReal& density,
    LemMatrix& primitives, LemMatrix& transport, LemMatrix& rhs,
    amrex::Real areaFactor, amrex::Real dt)
{
    TemperatureDiffusionScheme<DiffusionScheme::O2_Euler>(
        mass, density, primitives, transport, rhs, areaFactor, dt);
}

// -----------------------------------------------------------------------
// Public entry point: single-line dispatch to the selected specialisation
// -----------------------------------------------------------------------

void diffusion::TemperatureDiffusionImplementation(
    LemReal& mass, LemReal& density,
    LemMatrix& primitives, LemMatrix& transport, LemMatrix& rhs,
    const amrex::Real areaFactor, const amrex::Real dt)
{
    TemperatureDiffusionScheme<setup::diffusion_scheme>(
        mass, density, primitives, transport, rhs, areaFactor, dt);
}

// -----------------------------------------------------------------------
// Species diffusion (scheme-independent)
// -----------------------------------------------------------------------

void diffusion::SpeciesDiffusionImplementation(
    LemReal& mass, LemReal& density,
    LemMatrix& primitives, LemMatrix& transport, LemMatrix& rhs,
    const amrex::Real areaFactor)
{
    constexpr int N = NUM_LEM;
    auto eos = pele::physics::PhysicsType::eos();

    amrex::Real dm = 0.0;
    for (int i = 0; i < N; ++i) dm += mass[i];
    dm /= static_cast<amrex::Real>(N);

    std::array<std::array<amrex::Real, NUM_LEM>, NUM_SPECIES> X{};
    std::array<std::array<amrex::Real, NUM_LEM>, NUM_SPECIES> H{};

    for (int i = 0; i < N; ++i) {
        amrex::Real Y[NUM_SPECIES];
        for (int sp = 0; sp < NUM_SPECIES; ++sp) Y[sp] = primitives[i][sp];
        amrex::Real Xloc[NUM_SPECIES], Hloc[NUM_SPECIES];
        eos.Y2X(Y, Xloc);
        eos.T2Hi(primitives[i][NUM_SPECIES], Hloc);
        for (int sp = 0; sp < NUM_SPECIES; ++sp) { X[sp][i] = Xloc[sp]; H[sp][i] = Hloc[sp]; }
    }

    std::array<std::array<amrex::Real, NUM_LEM>, NUM_SPECIES> Ydot{};
    LemReal Tdot{}, fluxR{}, fluxL{};

    for (int i = 0; i < N; ++i) {
        const int iR = (i + 1) % N;
        const int iL = (i - 1 + N) % N;
        const amrex::Real rho  = density[i];
        const amrex::Real rhoR = density[iR];
        const amrex::Real rhoL = density[iL];

        for (int sp = 0; sp < NUM_SPECIES; ++sp) {
            const amrex::Real Xc = X[sp][i], XR = X[sp][iR], XL = X[sp][iL];
            const amrex::Real Yc = primitives[i][sp], YR = primitives[iR][sp], YL = primitives[iL][sp];

            const amrex::Real alphaR = SafeHarmonicMean(rhoR*transport[iR][sp], rho*transport[i][sp]);
            const amrex::Real alphaL = SafeHarmonicMean(rhoL*transport[iL][sp], rho*transport[i][sp]);

            const amrex::Real dX_R = (XR - Xc) / dm;
            const amrex::Real dX_L = (Xc - XL) / dm;

            fluxR[i] += alphaR * dX_R;
            fluxL[i] += alphaL * dX_L;

            Ydot[sp][i] = alphaR * dX_R - alphaL * dX_L;
            Tdot[i] += alphaR * dX_R * 0.25*(H[sp][iR]+H[sp][i])*(Yc+YR)
                     - alphaL * dX_L * 0.25*(H[sp][iL]+H[sp][i])*(Yc+YL);
        }
    }

    for (int i = 0; i < N; ++i) {
        const int iR = (i + 1) % N;
        const int iL = (i - 1 + N) % N;
        for (int sp = 0; sp < NUM_SPECIES; ++sp) {
            const amrex::Real corr = 0.5*(primitives[iR][sp]+primitives[i][sp])*fluxR[i]
                                   - 0.5*(primitives[iL][sp]+primitives[i][sp])*fluxL[i];
            Ydot[sp][i] = Ydot[sp][i] - corr;
        }
    }

    const amrex::Real scale = areaFactor * areaFactor / dm;
    for (int i = 0; i < N; ++i) {
        for (int sp = 0; sp < NUM_SPECIES; ++sp) rhs[i][sp] += Ydot[sp][i] * scale;
        rhs[i][NUM_SPECIES] += Tdot[i] * scale;
    }
}

// -----------------------------------------------------------------------
// Transport coefficients
// -----------------------------------------------------------------------

void diffusion::CalculationOfTransport(
    const TransParm* ltransparm,
    LemMatrix& transport, LemReal& density, LemMatrix& primitives)
{
    auto trans = pele::physics::PhysicsType::transport();
    auto eos   = pele::physics::PhysicsType::eos();

    for (int lem = 0; lem < NUM_LEM; lem++) {
        const amrex::Real rho = density[lem];
        const amrex::Real T   = primitives[lem][NUM_SPECIES];
        amrex::Real Y[NUM_SPECIES] = {0.0};
        for (int sp = 0; sp < NUM_SPECIES; sp++) Y[sp] = primitives[lem][sp];
        //Here D = rho*D
        amrex::Real Ddiag[NUM_SPECIES]={0.0}, chi[NUM_SPECIES]={0.0};
        amrex::Real mu, xi, lam;
        trans.transport(false, false, true, true, false, T, rho, Y, Ddiag, chi, mu, xi, lam, ltransparm);

        amrex::Real cp;
        eos.TY2Cp(T, Y, cp);

        for (int sp = 0; sp < NUM_SPECIES; sp++) transport[lem][sp] = Ddiag[sp];
        transport[lem][NUM_SPECIES] = lam / cp;   // rho*alpha_T = lambda/Cp
    }
}

// -----------------------------------------------------------------------
// Combined entry point — unified coupled sub-step for T and Y
//
// Structure:
//   1. CalculationOfTransport
//   2. Single dt_cfl = min stable step covering both thermal and species diffusivities
//   3. N_sub = ceil(dt / dt_cfl)  — shared by T and Y
//   4. Coupled Euler CDS loop:  advance T and Y together each sub-step
//   5. Write time-averaged rates to rhs
// -----------------------------------------------------------------------

void diffusion::ComputationOfTemperatureAndSpeciesDiffusion(
    LemReal& mass, LemReal& density,
    LemMatrix& primitives, LemMatrix& transport, LemMatrix& rhs,
    const amrex::Real areaFactor, const TransParm* ltransparm, const amrex::Real dt)
{
    constexpr int N = NUM_LEM;
    auto eos = pele::physics::PhysicsType::eos();

    // 1. Transport coefficients
    CalculationOfTransport(ltransparm, transport, density, primitives);

    // 2. dm, CFL timestep (ComputeStableTimestep already scans all diffusivities)
    amrex::Real dm = 0.0;
    for (int i = 0; i < N; i++) dm += mass[i];
    dm /= static_cast<amrex::Real>(N);

    const amrex::Real scale  = areaFactor * areaFactor / dm;
    const amrex::Real dt_cfl = ComputeStableTimestep(density, transport, dm, areaFactor, 0.5);
    const int         N_sub  = std::max(1, static_cast<int>(std::ceil(dt / dt_cfl)));
    const amrex::Real dt_sub = dt / static_cast<amrex::Real>(N_sub);

    // 3. Freeze X (mole fractions) and H (enthalpies) at T^n.
    //    X depends on Y; since Y changes slowly between sub-steps this is a good
    //    approximation and avoids repeated EOS calls inside the loop.
    std::array<std::array<amrex::Real, NUM_SPECIES>, NUM_LEM> X{}, H{};
    for (int i = 0; i < N; i++) {
        amrex::Real Y[NUM_SPECIES];
        for (int sp = 0; sp < NUM_SPECIES; sp++) Y[sp] = primitives[i][sp];
        amrex::Real Xloc[NUM_SPECIES], Hloc[NUM_SPECIES];
        eos.Y2X(Y, Xloc);
        eos.T2Hi(primitives[i][NUM_SPECIES], Hloc);
        for (int sp = 0; sp < NUM_SPECIES; sp++) { X[i][sp] = Xloc[sp]; H[i][sp] = Hloc[sp]; }
    }

    // 4. Scratch — working copy and accumulators
    LemMatrix& prim_work  = tl_scratch.prim_work;   prim_work  = primitives;
    LemReal&   src_T_sum  = tl_scratch.src_T_sum;   src_T_sum  = {};
    LemMatrix& src_Y_sum  = tl_scratch.src_Y_sum;   src_Y_sum  = {};
    LemReal&   src_T      = tl_scratch.src_T;
    LemMatrix& src_Y      = tl_scratch.src_Y;
    LemReal&   src_T_sp   = tl_scratch.src1;        // reuse slot for T contribution from species
    LemReal&   fluxR      = tl_scratch.fluxR;
    LemReal&   fluxL      = tl_scratch.fluxL;

    // 5. Coupled Euler CDS sub-step loop
    for (int sub = 0; sub < N_sub; sub++) {
        // Thermal diffusion source (O2 CDS on T)
        ComputeEnergySource_O2(density, prim_work, transport, dm, src_T);

        // Species diffusion sources (O2 CDS on Y, with enthalpy-flux T correction)
        ComputeSpeciesSource(density, prim_work, transport, dm, X, H,
                             src_Y, src_T_sp, fluxR, fluxL);

        // Advance prim_work: T gets both thermal and enthalpy-flux contributions
        for (int i = 0; i < N; i++) {
            prim_work[i][NUM_SPECIES] += dt_sub * (src_T[i] + src_T_sp[i]) * scale;
            for (int sp = 0; sp < NUM_SPECIES; sp++)
                prim_work[i][sp] += dt_sub * src_Y[i][sp] * scale;
        }

        // Accumulate time-averaged rates
        for (int i = 0; i < N; i++) {
            src_T_sum[i] += src_T[i] + src_T_sp[i];
            for (int sp = 0; sp < NUM_SPECIES; sp++)
                src_Y_sum[i][sp] += src_Y[i][sp];
        }
    }

    // 6. Write averaged rates to rhs
    const amrex::Real inv_N = 1.0 / static_cast<amrex::Real>(N_sub);
    for (int i = 0; i < N; i++) {
        rhs[i][NUM_SPECIES] += src_T_sum[i] * inv_N * scale;
        for (int sp = 0; sp < NUM_SPECIES; sp++)
            rhs[i][sp] += src_Y_sum[i][sp] * inv_N * scale;
    }
}

} // namespace clem
