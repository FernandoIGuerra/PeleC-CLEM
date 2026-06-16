#ifndef PARTICLE_FUNCTIONALITIES_H
#define PARTICLE_FUNCTIONALITIES_H

#include "ClemParticles/clem_particle.H"

namespace clem{
    namespace particle_functionalities{
        void SetParticleWithDefaultValues(  amrex::IntVect iv, const amrex::Real* plo, 
                                            const amrex::Real* dx, const amrex::Real dx_inner, 
                                            int index_lem, ClemParticles::ParticleType& particle );
        
        void SetParticleWithOtherParticleValuesAdvection(   amrex::IntVect& iv, const amrex::Real* plo, const amrex::Real* dx, 
                                                            const amrex::Real flux_out,
                                                            const amrex::Real net_flux,
                                                            const ClemParticles::ParticleType& p_in, 
                                                            ClemParticles::ParticleType& p_out);

        void CopyParticleDataToArray(std::array<amrex::Real, NUM_LEM>& mass,
                                    std::array<amrex::Real, NUM_LEM>& density,
                                    std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& primitives,
                                    std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& conservatives,
                                    std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& rhs,
                                    ClemParticles::ParticleType& p,
                                    int lem_index);
        
        void CopyArrayDataToParticle(const std::array<amrex::Real, NUM_LEM>& mass,
                                    const std::array<amrex::Real, NUM_LEM>& density,
                                    const std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& primitives,
                                    const std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& conservatives,
                                    const std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& rhs,
                                    ClemParticles::ParticleType& particle,
                                    int lem_index);
    }

}

#endif /* PARTICLE_FUNCTIONALITIES_H */