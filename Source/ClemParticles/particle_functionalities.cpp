#include "ClemParticles/particle_functionalities.H"
#include "ClemParticles/int_real_data.H"
#include "ClemParticles/clem_index_definition.H"


#include <AMReX.H>
#include <AMReX_ParticleContainer.H>
#include <AMReX_Particles.H>
#include <AMReX_ParIter.H>

namespace clem {
    void particle_functionalities::SetParticleWithDefaultValues(   amrex::IntVect iv, const amrex::Real* plo, 
                                                                    const amrex::Real* dx, const amrex::Real dx_inner, 
                                                                    int index_lem, ClemParticles::ParticleType& particle )
    {
        AMREX_D_TERM(
            particle.pos(0) = plo[0] + iv[0] * dx[0] + (index_lem + 1) * dx_inner;,
            particle.pos(1) = plo[1] + (iv[1] + 0.5) * dx[1];,
            particle.pos(2) = plo[2] + (iv[2] + 0.5) * dx[2];
        )

        AMREX_D_TERM(
            particle.idata(IntData::i) = iv[0];,
            particle.idata(IntData::j) = iv[1];,
            particle.idata(IntData::k) = iv[2];
        )   
        particle.idata(IntData::idx)      = index_lem;


        particle.idata(IntData::advected) = 0;

        //Definition of Real component 
        particle.rdata(RealData::mass)     = 0.0;
        particle.rdata(RealData::pressure) = 0.0;
        particle.rdata(RealData::volume)   = 0.0;
        particle.rdata(RealData::flux)     = 0.0;

        particle.rdata(RealData::density)  = 0.0;
        particle.rdata(RealData::energy)   = 0.0;
        particle.rdata(RealData::temp)     = 0.0;

        //Initialization of SP values
        for(unsigned int sp = 0; sp < NUM_SPECIES ; sp++){
            particle.rdata(RealData::SP(sp)) = 0.0;
        }
    }

    void particle_functionalities::SetParticleWithOtherParticleValuesAdvection(   amrex::IntVect& iv, const amrex::Real* plo, const amrex::Real* dx, 
                                                                                const amrex::Real flux_out,
                                                                                const amrex::Real net_flux,
                                                                                const ClemParticles::ParticleType& p_in, 
                                                                                ClemParticles::ParticleType& p_out)
    {
        AMREX_D_TERM(
            p_out.pos(0) = plo[0] + (iv[0] + 0.5) * dx[0] ;,
            p_out.pos(1) = plo[1] + (iv[1] + 0.5) * dx[1];,
            p_out.pos(2) = plo[2] + (iv[2] + 0.5) * dx[2];
        )

        p_out.idata(IntData::idx)      = 100; //Default value for the index of the lem, this is not important since will be updated later
        AMREX_D_TERM(
            p_out.idata(IntData::i) = iv[0];,
            p_out.idata(IntData::j) = iv[1];,
            p_out.idata(IntData::k) = iv[2];
        )   

        p_out.idata(IntData::advected) = 1;

        //Definition of Real component 
        p_out.rdata(RealData::mass)     = flux_out;
        p_out.rdata(RealData::pressure) = p_in.rdata(RealData::pressure);
        p_out.rdata(RealData::flux)     = std::abs(net_flux);

        p_out.rdata(RealData::density)  = p_in.rdata(RealData::density);
        p_out.rdata(RealData::energy)   = p_in.rdata(RealData::energy);
        p_out.rdata(RealData::temp)     = p_in.rdata(RealData::temp);

        //to preserve the density of the particle we set the volume based on the mass and density
        p_out.rdata(RealData::volume)   = p_out.rdata(RealData::mass) / p_out.rdata(RealData::density);
        //Initialization of SP values
        for(unsigned int sp = 0; sp < NUM_SPECIES ; sp++){
            p_out.rdata(RealData::SP(sp)) = p_in.rdata(RealData::SP(sp));
         }
    }



    void particle_functionalities::CopyParticleDataToArray( std::array<amrex::Real, NUM_LEM>& mass,
                                                            std::array<amrex::Real, NUM_LEM>& density,
                                                            std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& primitives,
                                                            std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& conservatives,
                                                            std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& rhs,
                                                            ClemParticles::ParticleType& particle,
                                                            int lem_index) 
    {
        const amrex::Real rho       = particle.rdata(RealData::density);
        const amrex::Real one_rho   = 1/rho;
        const amrex::Real m         = particle.rdata(RealData::mass);

        for( int sp = 0; sp < NUM_SPECIES; sp++){
            primitives[lem_index][sp]             = particle.rdata(RealData::SP(sp)) * one_rho;
            conservatives[lem_index][sp]          = particle.rdata(RealData::SP(sp));
        }

        mass[lem_index]                           = m;
        density[lem_index]                        = rho;
        primitives[lem_index][NUM_SPECIES]        = particle.rdata(RealData::temp);
        conservatives[lem_index][NUM_SPECIES]     = particle.rdata(RealData::energy);

    }

    void particle_functionalities::CopyArrayDataToParticle( const std::array<amrex::Real, NUM_LEM>& mass,
                                                            const std::array<amrex::Real, NUM_LEM>& density,
                                                            const std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& primitives,
                                                            const std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& conservatives,
                                                            const std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM>& rhs,
                                                            ClemParticles::ParticleType& particle,
                                                            int lem_index) 
    {

        amrex::Real rho = 0.0;
        for(int sp = 0; sp < NUM_SPECIES; sp++){
            particle.rdata(RealData::SP(sp))    = conservatives[lem_index][sp];
            rho                                 += conservatives[lem_index][sp];
        }
        particle.rdata(RealData::temp)      = primitives[lem_index][NUM_SPECIES];
        particle.rdata(RealData::energy)    = conservatives[lem_index][NUM_SPECIES];
        particle.rdata(RealData::density)   = rho;    
    }
}