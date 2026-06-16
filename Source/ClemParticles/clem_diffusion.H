#ifndef CLEM_DIFFUSION_H
#define CLEM_DIFFUSION_H

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>
#include <TransportParams.H>


#include "ClemParticles/clem_index_definition.H"
#include "ClemParticles/clem_setup.H"

#include <array>
//Implementation of diffussion
    /*
    @brief Implementation of Temeprature diffusion using centered finitie difference. all is accumualtes in rhs container
        rhs container has NUM_SPECIES + 1 values. The first NUM_SPECIES are for the mass fraction and the lastes is for the energy diffusion
        For both of diffusion, a periodic BC condition is set within the Clem array
    */

namespace clem{
    namespace diffusion{

        // Convenience aliases to reduce repetition in declarations
        using LemReal   = std::array<amrex::Real, NUM_LEM>;
        using LemMatrix = std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>;

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

        /**
         * Primary template for the thermal diffusion scheme.
         * Each DiffusionScheme value is fully specialised in clem_diffusion.cpp.
         * To add a new scheme: add an enum value in clem_numerics_config.H,
         * write a new explicit specialisation here, and select it in clem_setup.H.
         */
        template<DiffusionScheme Scheme>
        void TemperatureDiffusionScheme(
                                            LemReal&   mass,
                                            LemReal&   density,
                                            LemMatrix& primitives,
                                            LemMatrix& transport,
                                            LemMatrix& rhs,
                                            amrex::Real areaFactor,
                                            amrex::Real dt);

        // Public entry point — dispatches to TemperatureDiffusionScheme<setup::diffusion_scheme>
        void    TemperatureDiffusionImplementation(
                                            LemReal&   mass,
                                            LemReal&   density,
                                            LemMatrix& primitives,
                                            LemMatrix& transport,
                                            LemMatrix& rhs,
                                            const amrex::Real areaFactor,
                                            const amrex::Real dt = 0.0);


        void    SpeciesDiffusionImplementation(
                                            LemReal&   mass,
                                            LemReal&   density,
                                            LemMatrix& primitives,
                                            LemMatrix& transport,
                                            LemMatrix& rhs,
                                            const amrex::Real areaFactor);

        using TransParm = pele::physics::transport::TransParm<
                              pele::physics::PhysicsType::eos_type,
                              pele::physics::PhysicsType::transport_type>;

        void CalculationOfTransport(const TransParm* ltransparm,
                                    LemMatrix& transport,
                                    LemReal&   density,
                                    LemMatrix& primitives);

        void ComputationOfTemperatureAndSpeciesDiffusion(
                                            LemReal&   mass,
                                            LemReal&   density,
                                            LemMatrix& primitives,
                                            LemMatrix& transport,
                                            LemMatrix& rhs,
                                            const amrex::Real areaFactor,
                                            const TransParm* ltransparm,
                                            const amrex::Real dt = 0.0);
    }
}

#endif /* CLEM_DIFFUSION_H */