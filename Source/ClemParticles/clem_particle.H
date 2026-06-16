#ifndef CLEM_PARTICLES_H
#define CLEM_PARTICLES_H

#include <iostream>

#include <AMReX.H>
#include <AMReX_Particles.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParticleContainer.H>
#include <AMReX_BCRec.H>
#include <ReactorBase.H>
#include "PelePhysics.H"
#include "mechanism.H"

#include "ClemParticles/clem_index_definition.H"
#include "ClemParticles/int_real_data.H"

namespace clem{

//==================================================
// Particle container
//==================================================

    class ClemParticles    
        : public amrex::ParticleContainer<RealData::ncomps, IntData::ncomps>
    {
        public:
            using ParticleType = typename amrex::ParticleContainer<RealData::ncomps, IntData::ncomps>::ParticleType;
            using MyPartIter = amrex::ParIter<RealData::ncomps, IntData::ncomps>;
            //Initializatior
            ClemParticles();
            ClemParticles(  const amrex::Geometry& geom, 
                            const amrex::DistributionMapping& dm, 
                            const amrex::BoxArray& grid);
            
            //Particle Generator
            void InitParticles();

            void DefineOwnershipParticleMesh();
            
            void UpdateIndexAfterAdvection();

            void SetParticlePropertiesFromMF(const amrex::MultiFab& mf);

            void GetCurrentAmountOfMassSpecies(amrex::MultiFab& mf);
            
            //inflow implementation - state_bc is the FillPatched LES state whose ghost cells (bcnormal) supply the incoming-fluid properties
            void ParticlesAdvection(amrex::Real dt, const amrex::MultiFab& state_bc);

            //outflow implementation - store the physical BC types (PCPhysBCType encoding) so advection knows where mass leaves the domain
            void SetPhysBC(const amrex::BCRec& bc) {
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    m_lo_bc[dir] = bc.lo(dir);
                    m_hi_bc[dir] = bc.hi(dir);
                }
            }

            void MoveParticles();
            
            void DoIsentropicPressureChange(amrex::MultiFab& Pressure);
            
            ParticleType CreateParticleTypeDefault(int index_lem, amrex::IntVect iv, const amrex::Real* plo, const amrex::Real* dx, amrex::Real dx_inner);

            ParticleType CreateParticleAdvectionTypeDefault(amrex::IntVect iv, const amrex::Real* plo, const amrex::Real* dx, const amrex::Real flux_out, const ParticleType& p_in);

            
            void CreateParticleBeforeRegridding();
               
            void DoRegridding();

            void DoClem(pele::physics::transport::TransParm<pele::physics::PhysicsType::eos_type,pele::physics::PhysicsType::transport_type> const* ltransparm,
                        pele::physics::reactions::ReactorBase* reactor,
                        amrex::MultiFab& react_src,
                        amrex::Real dt,
                        amrex::MultiFab& nonReactSrc,
                        amrex::MultiFab& S_old,
                        amrex::MultiFab& S_new,
                        //inflow implementation - FillPatched state (ghost cells hold the bcnormal inflow Dirichlet values)
                        const amrex::MultiFab& Sborder);
            
            amrex::MultiFab& getFluxContainer(){return FluxFaces;}
            
        private:
            bool m_cell_vector_updated = false;
            // For each grid: BaseFab indexed by cell, holding particle indices
            std::map<int, amrex::BaseFab<std::vector<int>>> m_cell_vectors;
            amrex::MultiFab FluxFaces;
            //outflow implementation - physical BC type per direction (PCPhysBCType encoding); default 0 = interior
            amrex::IntVect m_lo_bc{AMREX_D_DECL(0, 0, 0)};
            amrex::IntVect m_hi_bc{AMREX_D_DECL(0, 0, 0)};
    };


}//namespace

#endif //definition