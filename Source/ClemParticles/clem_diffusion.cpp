
#include "ClemParticles/clem_diffusion.H"
#include <TransportParams.H>
#include "ClemParticles/clem_index_definition.H"

namespace clem{
        void diffusion::TemperatureDiffusionImplementation(
                                            std::array<amrex::Real, NUM_LEM>& mass,
                                            std::array<amrex::Real, NUM_LEM>& density,
                                            std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& primitives,
                                            std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& transport,
                                            std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& rhs,
                                            const amrex::Real areaFactor)
    {
        constexpr int N = NUM_LEM;
        // Compute average mass factor (like dm)
        amrex::Real dm = 0.0;
        for (int i = 0; i < N; i++) dm  += mass[i];
        dm                              /= static_cast<amrex::Real>(N);

        std::vector<amrex::Real> EintSource(N, 0.0);

        // Loop over LEMs and compute thermal diffusion
        for (int i = 0; i < N; i++) {
            int iR = (i + 1) % N;       // right neighbor
            int iL = (i - 1 + N) % N;   // left neighbor

            // Densities
            amrex::Real rho  = density[i];
            amrex::Real rhoR = density[iR];
            amrex::Real rhoL = density[iL];

            // Temperatures
            amrex::Real T  = primitives[i][NUM_SPECIES];
            amrex::Real TR = primitives[iR][NUM_SPECIES];
            amrex::Real TL = primitives[iL][NUM_SPECIES];

            // Transport coefficients
            amrex::Real gamma  = rho*transport[i][CLEM_dComp_lambda];
            amrex::Real gammaR = rhoR*transport[iR][CLEM_dComp_lambda];
            amrex::Real gammaL = rhoL*transport[iL][CLEM_dComp_lambda];

            // Harmonic-averaged diffusivity
            amrex::Real alphaR = SafeHarmonicMean(gammaR, gamma);
            amrex::Real alphaL = SafeHarmonicMean(gammaL, gamma);

            // Fluxes
            amrex::Real F_R = alphaR * (TR - T) / dm;
            amrex::Real F_L = alphaL * (T - TL) / dm;

            // Store energy source
            EintSource[i] = rho * (F_R - F_L);
        }

        for (int i = 0; i < N; i++) {
            rhs[i][NUM_SPECIES] += EintSource[i] * areaFactor * areaFactor / dm;
        }  
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
                                            pele::physics::transport::TransParm<pele::physics::PhysicsType::eos_type,pele::physics::PhysicsType::transport_type> const* ltransparm)
    {
        CalculationOfTransport(ltransparm, transport, density, primitives);
        TemperatureDiffusionImplementation(mass, density, primitives, transport, rhs, areaFactor);
        SpeciesDiffusionImplementation(mass,density,primitives,transport,rhs,areaFactor);
    }
    
}