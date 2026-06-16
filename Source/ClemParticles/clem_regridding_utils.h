#ifndef CLEM_REGRIDDING_UTILS_H
#define CLEM_REGRIDDING_UTILS_H

#include <AMReX_REAL.H>
#include <AMReX.H>
#include <array>

namespace clem {
        /*
        -> internal energy (ueint) is computed one the mass average of temperature is done ueint(T) this is done outside this function
        -> after the regridding, we have some particles with mass 0.0, those must be with index higher or equal to NUM_LEM
        -> p.rdata(RealData::SP(sp)) contains conservative form rho*y
        ->Pressure is mass weighted intentionally. Note that P is set equal in all elements previously, so no change.

        */
    namespace regridding {

        struct LemAccumulator
        {
            amrex::Real mass      = 0.0;
            amrex::Real pressure  = 0.0;
            amrex::Real volume    = 0.0;
            amrex::Real density   = 0.0;
            amrex::Real energy    = 0.0;
            amrex::Real temp      = 0.0;
            std::array<amrex::Real, NUM_SPECIES> sp{0.0};

            AMREX_GPU_HOST_DEVICE
            void reset();
        };

        // --- Function declarations (inline for GPU safety)

        AMREX_GPU_HOST_DEVICE inline
        void AccumulateDataIntoStructureGivenMassPortion(
                                                            LemAccumulator& acc,
                                                            const ClemParticles::ParticleType& p,
                                                            amrex::Real mass_need);

        AMREX_GPU_HOST_DEVICE inline
        void AccumulateDataIntoStructureGivenMassPortion(
                                                            LemAccumulator& acc,
                                                            const LemAccumulator& aux,
                                                            amrex::Real mass_need);

        AMREX_GPU_HOST_DEVICE inline
        void StoreRemainningDataIntoAux(
                                        LemAccumulator& data,
                                        const ClemParticles::ParticleType& p);

        AMREX_GPU_HOST_DEVICE inline
        void UpdateLemFromAccumulator(
                                        const LemAccumulator& data,
                                        ClemParticles::ParticleType& dst,
                                        amrex::Real new_dm);

        AMREX_GPU_HOST_DEVICE inline
        void AccumulateDataIntoStructure(
                                        LemAccumulator& acc,
                                        LemAccumulator& aux);

        AMREX_GPU_HOST_DEVICE inline
        void AccumulateDataIntoStructure(
                                        LemAccumulator& acc,
                                        ClemParticles::ParticleType& p);

    } // namespace RegriddingUt

}//namespace clem
#endif