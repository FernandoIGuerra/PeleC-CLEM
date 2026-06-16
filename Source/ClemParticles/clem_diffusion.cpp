
#include "ClemParticles/clem_diffusion.H"
#include <TransportParams.H>
#include "ClemParticles/clem_index_definition.H"
#include <limits>
#include <cmath>

namespace clem{

// -----------------------------------------------------------------------
// Internal helpers — one per spatial/time scheme.  Only the one selected
// by CLEM_DIFF_SCHEME is compiled into the final binary.
// -----------------------------------------------------------------------

#if CLEM_DIFF_SCHEME == CLEM_DIFF_SCHEME_O2_EULER   || \
    CLEM_DIFF_SCHEME == CLEM_DIFF_SCHEME_O2_RK2     || \
    CLEM_DIFF_SCHEME == CLEM_DIFF_SCHEME_O2_SUBSTEP

// Compute the per-particle energy flux divergence using the 2nd-order
// central-difference stencil.  Result written into EintSource[N].
static void ComputeEnergySource_O2(
    const std::array<amrex::Real, NUM_LEM>&                                density,
    const std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>&   T_field,
    const std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>&   transport,
    amrex::Real dm,
    std::array<amrex::Real, NUM_LEM>&                                       EintSource)
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

        const amrex::Real gamma  = rho  * transport[i][CLEM_dComp_lambda];
        const amrex::Real gammaR = rhoR * transport[iR][CLEM_dComp_lambda];
        const amrex::Real gammaL = rhoL * transport[iL][CLEM_dComp_lambda];

        const amrex::Real alphaR = diffusion::SafeHarmonicMean(gammaR, gamma);
        const amrex::Real alphaL = diffusion::SafeHarmonicMean(gammaL, gamma);

        // 2nd-order face fluxes:  F = alpha * dT/dm
        const amrex::Real F_R = alphaR * (TR - T) / dm;
        const amrex::Real F_L = alphaL * (T - TL) / dm;

        EintSource[i] = rho * (F_R - F_L);
    }
}
#endif // O2 stencil needed

#if CLEM_DIFF_SCHEME == CLEM_DIFF_SCHEME_O4_EULER

// Compute the per-particle energy flux divergence using the 4th-order
// accurate face-flux stencil (conservative form).
//
// 4th-order first derivative at face i+1/2:
//   (dT/dm) = (-T_{i+2} + 27T_{i+1} - 27T_i + T_{i-1}) / (24*dm)
//
// Stability note: the maximum eigenvalue of the O4 Laplacian is larger
// than the O2 one, so the stability limit tightens to r <= 3/11 ≈ 0.27.
static void ComputeEnergySource_O4(
    const std::array<amrex::Real, NUM_LEM>&                                density,
    const std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>&   T_field,
    const std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>&   transport,
    amrex::Real dm,
    std::array<amrex::Real, NUM_LEM>&                                       EintSource)
{
    constexpr int N = NUM_LEM;
    constexpr amrex::Real c24 = 1.0 / 24.0;

    for (int i = 0; i < N; i++) {
        const int iR  = (i + 1) % N;
        const int iRR = (i + 2) % N;
        const int iL  = (i - 1 + N) % N;
        const int iLL = (i - 2 + N) % N;

        const amrex::Real rho  = density[i];
        const amrex::Real rhoR = density[iR];
        const amrex::Real rhoL = density[iL];

        const amrex::Real T   = T_field[i][NUM_SPECIES];
        const amrex::Real TR  = T_field[iR][NUM_SPECIES];
        const amrex::Real TRR = T_field[iRR][NUM_SPECIES];
        const amrex::Real TL  = T_field[iL][NUM_SPECIES];
        const amrex::Real TLL = T_field[iLL][NUM_SPECIES];

        // Harmonic-mean transport at nearest-neighbour faces (same as O2)
        const amrex::Real gamma  = rho  * transport[i][CLEM_dComp_lambda];
        const amrex::Real gammaR = rhoR * transport[iR][CLEM_dComp_lambda];
        const amrex::Real gammaL = rhoL * transport[iL][CLEM_dComp_lambda];

        const amrex::Real alphaR = diffusion::SafeHarmonicMean(gammaR, gamma);
        const amrex::Real alphaL = diffusion::SafeHarmonicMean(gammaL, gamma);

        // 4th-order face gradient:  (-T_{i+2}+27T_{i+1}-27T_i+T_{i-1})/(24*dm)
        const amrex::Real F_R = alphaR * (-TRR + 27.0*TR - 27.0*T + TL)  * c24 / dm;
        const amrex::Real F_L = alphaL * (-TR  + 27.0*T  - 27.0*TL + TLL)* c24 / dm;

        EintSource[i] = rho * (F_R - F_L);
    }
}
#endif // O4 stencil needed

#if CLEM_DIFF_SCHEME == CLEM_DIFF_SCHEME_O2_SUBSTEP

// Compute the CFL-stable diffusion sub-step size.
//
// beta = max_i( rho_i * max(lambda/Cp, D_sp) ) — the largest effective
// mass-space diffusivity weighted by density, matching the stability analysis
// of the 3-point explicit stencil:  dt_cfl = dm^2 / (2 * beta).
//
// transport[i][sp]            = D_sp      (or rho*D_sp from PelePhysics)
// transport[i][CLEM_dComp_lambda] = lambda/Cp = rho*alpha_T
// Both already carry the rho factor so we multiply by density[i] again
// to get rho^2*alpha, matching the "beta = D*rho^2" formula in getTimeCFL.
static amrex::Real ComputeStableTimestep(
    const std::array<amrex::Real, NUM_LEM>&                                density,
    const std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>&   transport,
    amrex::Real dm)
{
    amrex::Real beta = 0.0;
    for (int i = 0; i < NUM_LEM; i++) {
        amrex::Real alpha_max = transport[i][CLEM_dComp_lambda];
        for (int sp = 0; sp < NUM_SPECIES; sp++)
            alpha_max = std::max(alpha_max, transport[i][sp]);
        beta = std::max(beta, alpha_max * density[i]);
    }
    return (beta > 0.0) ? (dm * dm) / (2.0 * beta)
                        : std::numeric_limits<amrex::Real>::max();
}

#endif // SUBSTEP helper

// -----------------------------------------------------------------------
// Public function: dispatch to the compile-time selected scheme
// -----------------------------------------------------------------------
    void diffusion::TemperatureDiffusionImplementation(
                                            std::array<amrex::Real, NUM_LEM>& mass,
                                            std::array<amrex::Real, NUM_LEM>& density,
                                            std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& primitives,
                                            std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& transport,
                                            std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& rhs,
                                            const amrex::Real areaFactor,
                                            const amrex::Real dt)
    {
        constexpr int N = NUM_LEM;

        amrex::Real dm = 0.0;
        for (int i = 0; i < N; i++) dm += mass[i];
        dm /= static_cast<amrex::Real>(N);

        const amrex::Real scale = areaFactor * areaFactor / dm;

        std::array<amrex::Real, NUM_LEM> EintSource{};

#if CLEM_DIFF_SCHEME == CLEM_DIFF_SCHEME_O2_EULER

        ComputeEnergySource_O2(density, primitives, transport, dm, EintSource);
        for (int i = 0; i < N; i++)
            rhs[i][NUM_SPECIES] += EintSource[i] * scale;

#elif CLEM_DIFF_SCHEME == CLEM_DIFF_SCHEME_O4_EULER

        ComputeEnergySource_O4(density, primitives, transport, dm, EintSource);
        for (int i = 0; i < N; i++)
            rhs[i][NUM_SPECIES] += EintSource[i] * scale;

#elif CLEM_DIFF_SCHEME == CLEM_DIFF_SCHEME_O2_RK2

        // Heun's method (explicit RK2) for the diffusion operator only.
        // Transport coefficients are frozen at T^n throughout both stages.
        //
        // Stage 1: stencil at T^n → EintSource1
        // Predictor: T^* = T^n + dt * EintSource1 * scale / (rho * Cv)
        // Stage 2: stencil at T^* → EintSource2
        // Corrector: rhs += 0.5*(EintSource1 + EintSource2) * scale

        auto eos = pele::physics::PhysicsType::eos();

        // --- Stage 1 ---
        std::array<amrex::Real, NUM_LEM> EintSource2{};
        ComputeEnergySource_O2(density, primitives, transport, dm, EintSource);

        // Build predictor temperature array (working copy of primitives)
        std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM> prim_star = primitives;

        for (int i = 0; i < N; i++) {
            amrex::Real massfrac[NUM_SPECIES];
            for (int sp = 0; sp < NUM_SPECIES; sp++)
                massfrac[sp] = primitives[i][sp];

            amrex::Real Cv = 0.0;
            eos.TY2Cv(primitives[i][NUM_SPECIES], massfrac, Cv);

            // dT/dt from the diffusion operator
            const amrex::Real dTdt = (Cv > 0.0)
                ? EintSource[i] * scale / (density[i] * Cv)
                : 0.0;

            prim_star[i][NUM_SPECIES] = primitives[i][NUM_SPECIES] + dt * dTdt;
        }

        // --- Stage 2 ---
        ComputeEnergySource_O2(density, prim_star, transport, dm, EintSource2);

        // --- Corrector: average the two stages ---
        for (int i = 0; i < N; i++)
            rhs[i][NUM_SPECIES] += 0.5 * (EintSource[i] + EintSource2[i]) * scale;

#elif CLEM_DIFF_SCHEME == CLEM_DIFF_SCHEME_O2_SUBSTEP

        // CFL-based sub-stepping with 2nd-order spatial stencil.
        //
        // Algorithm:
        //   1. Compute dt_cfl from the diffusion stability condition.
        //   2. Take N_sub = ceil(dt / dt_cfl) explicit Euler steps of size
        //      dt_sub = dt / N_sub using a LOCAL copy of the temperature field.
        //      Species and density are frozen (transport coefficients frozen too).
        //   3. At each sub-step accumulate EintSource.
        //   4. The averaged EintSource is passed to rhs so the reactor applies
        //      the correct total energy change over the full dt:
        //        delta(rho*e)_diff = dt * avg_rate = dt_sub * sum(EintSource_k) * scale
        //
        // This scheme is unconditionally stable for any dt and more accurate
        // than the single frozen-rate approach for large Courant numbers.

        auto eos = pele::physics::PhysicsType::eos();

        const amrex::Real dt_cfl = ComputeStableTimestep(density, transport, dm);
        const int N_sub = std::max(1, static_cast<int>(std::ceil(dt / dt_cfl)));
        const amrex::Real dt_sub = dt / static_cast<amrex::Real>(N_sub);

        // Working copy of primitives — only temperature column is updated
        std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM> prim_work = primitives;

        // Accumulated EintSource across all sub-steps
        std::array<amrex::Real, NUM_LEM> EintSource_sum{};

        for (int sub = 0; sub < N_sub; sub++) {
            std::array<amrex::Real, NUM_LEM> EintSource_sub{};
            ComputeEnergySource_O2(density, prim_work, transport, dm, EintSource_sub);

            for (int i = 0; i < N; i++) {
                // Advance local temperature: dT/dt = EintSource * scale / (rho * Cv)
                amrex::Real massfrac[NUM_SPECIES];
                for (int sp = 0; sp < NUM_SPECIES; sp++)
                    massfrac[sp] = primitives[i][sp];   // species frozen at T^n

                amrex::Real Cv = 0.0;
                eos.TY2Cv(prim_work[i][NUM_SPECIES], massfrac, Cv);

                const amrex::Real dTdt = (density[i] * Cv > 0.0)
                    ? EintSource_sub[i] * scale / (density[i] * Cv)
                    : 0.0;

                prim_work[i][NUM_SPECIES] += dt_sub * dTdt;
                EintSource_sum[i]         += EintSource_sub[i];
            }
        }

        // Pass the time-averaged rate so the reactor applies the correct delta(rho*e)
        const amrex::Real inv_Nsub = 1.0 / static_cast<amrex::Real>(N_sub);
        for (int i = 0; i < N; i++)
            rhs[i][NUM_SPECIES] += EintSource_sum[i] * inv_Nsub * scale;

#else
#error "Unknown CLEM_DIFF_SCHEME value — set 0, 1, 2, or 3 in clem_numerics_config.H"
#endif
    }

    void diffusion::SpeciesDiffusionImplementation(
                                            std::array<amrex::Real, NUM_LEM>& mass,
                                            std::array<amrex::Real, NUM_LEM>& density,
                                            std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& primitives,
                                            std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& transport,
                                            std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& rhs,
                                            const amrex::Real areaFactor)
    {
        constexpr int N = NUM_LEM;
        auto eos = pele::physics::PhysicsType::eos();

        amrex::Real dm = 0.0;
        for (int i = 0; i < N; ++i) dm += mass[i];
        dm /= static_cast<amrex::Real>(N);
        

        // Mole fractions and enthalpies
        std::array<std::array<amrex::Real, NUM_LEM>, NUM_SPECIES> X{};
        std::array<std::array<amrex::Real, NUM_LEM>, NUM_SPECIES> H{};

        for (int i = 0; i < N; ++i) {
            amrex::Real rho         = density[i];
            amrex::Real one_rho     = 1.0 / rho;

            amrex::Real Y[NUM_SPECIES];
            for (int sp = 0; sp < NUM_SPECIES; ++sp)
                Y[sp]               = primitives[i][sp];

            amrex::Real Xloc[NUM_SPECIES];
            eos.Y2X(Y, Xloc);
            amrex::Real Hloc[NUM_SPECIES];
            eos.T2Hi(primitives[i][NUM_SPECIES], Hloc);

            for (int sp = 0; sp < NUM_SPECIES; ++sp) {
                X[sp][i] = Xloc[sp];
                H[sp][i] = Hloc[sp];
            }
        }

        std::array<std::array<amrex::Real, NUM_LEM>, NUM_SPECIES> Ydot{};
        std::array<amrex::Real, NUM_LEM> Tdot{};
        std::array<amrex::Real, NUM_LEM> fluxR{}, fluxL{};

        // --- raw fluxes
        for (int i = 0; i < N; ++i) {
            int iR = (i + 1) % N;
            int iL = (i - 1 + N) % N;

            amrex::Real rho  = density[i];
            amrex::Real rhoR = density[iR];
            amrex::Real rhoL = density[iL];

            for (int sp = 0; sp < NUM_SPECIES; ++sp) {

                amrex::Real Xc  = X[sp][i];
                amrex::Real XR  = X[sp][iR];
                amrex::Real XL  = X[sp][iL];

                amrex::Real Yc  = primitives[i][sp];
                amrex::Real YR  = primitives[iR][sp];
                amrex::Real YL  = primitives[iL][sp];

                amrex::Real a  = rho  * transport[i][sp];
                amrex::Real aR = rhoR * transport[iR][sp];
                amrex::Real aL = rhoL * transport[iL][sp];

                amrex::Real alphaR = SafeHarmonicMean(aR, a);
                amrex::Real alphaL = SafeHarmonicMean(aL, a);

                amrex::Real dX_R = (XR - Xc) / dm;
                amrex::Real dX_L = (Xc - XL) / dm;

                amrex::Real F_R = alphaR * dX_R;
                amrex::Real F_L = alphaL * dX_L;

                fluxR[i] += F_R;
                fluxL[i] += F_L;

                Ydot[sp][i] = F_R - F_L;
                Tdot[i] += alphaR * dX_R * 0.25 * (H[sp][iR] + H[sp][i]) * (Yc + YR) - alphaL * dX_L * 0.25 * (H[sp][iL] + H[sp][i]) * (Yc + YL);
            }
        }

        // --- correction + conservative form
        for (int i = 0; i < N; ++i) {
            int iR = (i + 1) % N;
            int iL = (i - 1 + N) % N;

            amrex::Real rho = density[i];

            for (int sp = 0; sp < NUM_SPECIES; ++sp) {

                amrex::Real YR_ = 0.5 * (primitives[iR][sp] + primitives[i][sp]);
                amrex::Real YL_ = 0.5 * (primitives[iL][sp] + primitives[i][sp]);

                amrex::Real corr = YR_ * fluxR[i] - YL_ * fluxL[i];

                Ydot[sp][i] = (Ydot[sp][i] - corr) * rho;
            }

            Tdot[i] *= rho;
        }

        // --- RHS accumulation
        for (int i = 0; i < N; ++i) {
            for (int sp = 0; sp < NUM_SPECIES; ++sp) {
                rhs[i][sp] += Ydot[sp][i] * areaFactor * areaFactor / dm;
            }
            rhs[i][NUM_SPECIES] += Tdot[i] * areaFactor * areaFactor / dm;
        }
    }//function

    
    void diffusion::CalculationOfTransport(pele::physics::transport::TransParm<pele::physics::PhysicsType::eos_type,pele::physics::PhysicsType::transport_type> const* ltransparm,
                                            std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& transport,
                                            std::array<amrex::Real, NUM_LEM>& density,
                                            std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& primitives)
    {
        auto trans                  = pele::physics::PhysicsType::transport();
        auto eos                    = pele::physics::PhysicsType::eos();

        for(int lem = 0; lem < NUM_LEM; lem++){
            amrex::Real rho                         = density[lem];
            amrex::Real T                           = primitives[lem][NUM_SPECIES];
            amrex::Real massloc[NUM_SPECIES]       = {0.0};

            for(int sp = 0; sp < NUM_SPECIES ; sp++)
                massloc[sp] = primitives[lem][sp];
            
            //Data Conatiner to compute the transport variables
            amrex::Real Ddiag[NUM_SPECIES] = {0.0};
            amrex::Real chi_mix[NUM_SPECIES] = {0.0};
            amrex::Real muloc, xiloc, lamloc;

            //True only lambda(heat tranpor), and Ddiag (species diffusion)
            trans.transport(false, false, true, true, false,T, rho, massloc, Ddiag, chi_mix, muloc, xiloc, lamloc, ltransparm);

            amrex::Real cp;
            eos.TY2Cp(T, massloc, cp);
            //Storing data into 
            for(int sp = 0; sp < NUM_SPECIES; sp++)
                transport[lem][sp]  = Ddiag[sp];
            //Heat diffusivity [m**2/s]
            // heat diffusivity * rho // alpha = lambda /(cp * rho)
            //To be consitant with species givena s rho*D then we store also rho*alpha = lambda/cp
            transport[lem][NUM_SPECIES]     = lamloc/cp;
        }
    }

    void diffusion::ComputationOfTemperatureAndSpeciesDiffusion(
                                            std::array<amrex::Real, NUM_LEM>& mass,
                                            std::array<amrex::Real, NUM_LEM>& density,
                                            std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& primitives,
                                            std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& transport,
                                            std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& rhs,
                                            const amrex::Real areaFactor,
                                            pele::physics::transport::TransParm<pele::physics::PhysicsType::eos_type,pele::physics::PhysicsType::transport_type> const* ltransparm,
                                            const amrex::Real dt)
    {
        CalculationOfTransport(ltransparm, transport, density, primitives);
        //This is the dt super grid(LES)
        //Iteration over the CFL conditioned time usign teh subgrid properties.
        TemperatureDiffusionImplementation(mass, density, primitives, transport, rhs, areaFactor, dt);
        SpeciesDiffusionImplementation(mass, density, primitives, transport, rhs, areaFactor);
    }
    
}