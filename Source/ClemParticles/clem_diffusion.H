#ifndef CLEM_DIFFUSION_H
#define CLEM_DIFFUSION_H

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>
#include <TransportParams.H>


#include "ClemParticles/clem_index_definition.H"


#include <array>
//Implementation of diffussion
    /*
    @brief Implementation of Temeprature diffusion using centered finitie difference. all is accumualtes in rhs container
        rhs container has NUM_SPECIES + 1 values. The first NUM_SPECIES are for the mass fraction and the lastes is for the energy diffusion
        For both of diffusion, a periodic BC condition is set within the Clem array
    */

namespace clem{
    namespace diffusion{
        /**
         * @brief Compute the harmonic mean of two values, with a safeguard against very small values to prevent numerical instability.
         */
        AMREX_FORCE_INLINE
        amrex::Real SafeHarmonicMean(amrex::Real a, amrex::Real b){
            const amrex::Real eps = 1e-15;
            if (a * b > eps) {
                return 2.0 * a * b / (a + b);
            } else {
                return 0.0;
            }
        }

        void    TemperatureDiffusionImplementation(
                                            std::array<amrex::Real, NUM_LEM>& mass,
                                            std::array<amrex::Real, NUM_LEM>& density,
                                            std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& primitives,
                                            std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& transport,
                                            std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& rhs,
                                            const amrex::Real areaFactor);


        void    SpeciesDiffusionImplementation(
                                            std::array<amrex::Real, NUM_LEM>& mass,
                                            std::array<amrex::Real, NUM_LEM>& density,
                                            std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& primitives,
                                            std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& transport,
                                            std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& rhs,
                                            const amrex::Real areaFactor);
        
        void CalculationOfTransport(pele::physics::transport::TransParm<pele::physics::PhysicsType::eos_type,pele::physics::PhysicsType::transport_type> const* ltransparm,
                                    std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& transport,
                                    std::array<amrex::Real, NUM_LEM>& density,
                                    std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& primitives);

        /**
         * @brief High-level function to compute both temperature and species diffusion. It first calculates the transport properties, 
         * then applies the temperature diffusion and finally the species diffusion. The results are accumulated in the rhs container, 
         * which has NUM_SPECIES + 1 entries (the first NUM_SPECIES for mass fractions and the last one for energy diffusion).
         */

        void ComputationOfTemperatureAndSpeciesDiffusion(
                                            std::array<amrex::Real, NUM_LEM>& mass,
                                            std::array<amrex::Real, NUM_LEM>& density,
                                            std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& primitives,
                                            std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& transport,
                                            std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& rhs,
                                            const amrex::Real areaFactor,
                                            pele::physics::transport::TransParm<pele::physics::PhysicsType::eos_type,pele::physics::PhysicsType::transport_type> const* ltransparm);
    }
}

#endif /* CLEM_DIFFUSION_H */