#include "PeleC.H"
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>
#include <array>

#include <ReactorBase.H>

#include "ClemParticles/clem_particle.H"
#include "ClemParticles/int_real_data.H"
#include "ClemParticles/clem_index_definition.H"
#include "ClemParticles/clem_functionalities.H"
#include "ClemParticles/clem_advection.H"
#include "ClemParticles/particle_functionalities.H"
#include "ClemParticles/clem_diffusion.H"
//This libary is not used currently but it will be used in the future for the regridding process to avoid strcuture allcoation
#include "ClemParticles/clem_regridding_utils.h"

std::unique_ptr<clem::ClemParticles> PeleC::ClemContainer = nullptr;

namespace clem {
    ClemParticles::ClemParticles() {}
    //==================================================
    ClemParticles::ClemParticles  (const amrex::Geometry& geom,
                                        const amrex::DistributionMapping& dm,
                                        const amrex::BoxArray& grid)
        : amrex::ParticleContainer<RealData::ncomps, IntData::ncomps>(geom, dm, grid)
    {   
        constexpr int number_faces   = (AMREX_SPACEDIM == 2) ? 4 : 6; 
        FluxFaces.define(grid, dm, number_faces, 4);
    }

    
    /**
    @brief Particles allocation in the domain
    */
    void
    ClemParticles::InitParticles()
    {
        const int lev                   = 0;
        const auto& geom                = Geom(lev);
        const auto* dx                  = geom.CellSize();
        const auto* plo                 = geom.ProbLo();
        constexpr int numberParticles   = NUM_LEM;

        const amrex::Real dx_inner       = dx[0] / static_cast<amrex::Real>(numberParticles + 1);

        for (amrex::MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx    = mfi.tilebox();
            const int grid_id       = mfi.index();
            const int tile_id       = mfi.LocalTileIndex();
            auto& ptile             = GetParticles(lev)[{grid_id, tile_id}];
        
            for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv))
            {
                for (int ip = 0; ip < numberParticles; ++ip)
                {
                    ParticleType p;
                    p.id()      = ParticleType::NextID();
                    p.cpu()     = amrex::ParallelDescriptor::MyProc();
                    particle_functionalities::SetParticleWithDefaultValues(iv, plo, dx, dx_inner, ip, p);
                    ptile.push_back(p); 
            }
        }
        m_cell_vector_updated = false;
        }
    }
    /*
    @brief : PArticles are linked to the cell that they belong to
    */
    void
    ClemParticles::DefineOwnershipParticleMesh()
    {
        const int lev = 0;
        if (m_cell_vector_updated) return;

        for(amrex::MFIter mfi = MakeMFIter(lev, false); mfi.isValid(); ++mfi){
            const int grid_id           = mfi.index();
            auto& fab = m_cell_vectors[grid_id];
            fab.resize(mfi.validbox());
            fab.setVal(std::vector<int>{});
        }
        amrex::Print() << "Particle ownership mesh defined" << std::endl;
        //Here the particles are associated to a cell 
        for(MyPartIter pti(*this, lev); pti.isValid(); ++pti){
            const int np                    = pti.numParticles();
            const int grid_id               = pti.index();
            const int tile_id               = pti.LocalTileIndex();
            auto& particle_tile             = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
            auto& particles                 = particle_tile.GetArrayOfStructs();

            for(int pindex = 0; pindex < np; ++pindex)
            {
                //To a bf(i,j,k) is append the index of a particle
                ParticleType& p             = particles[pindex];
                if(p.rdata(RealData::mass) > 0.0001){
                    amrex::Print() << "Particle mass: " << p.rdata(RealData::mass) << " Particle ID: " << p.idata(IntData::idx) << " Position: " << p.idata(IntData::i) << ", " << p.idata(IntData::j) << std::endl;
                } 
                amrex::IntVect iv           = this->Index(p, lev); //Get the index from the position
                AMREX_ASSERT(m_cell_vectors[grid_id].box().contains(iv));
                m_cell_vectors[grid_id](iv).push_back(pindex); //add to each cell the belonging cells
            }

        }
        m_cell_vector_updated = true;
        //WritePlotFile("plt_particle_ownership", "particles");
        //This is a general algorthmin that sort the particles according to their index and flux 
        UpdateIndexAfterAdvection();
    }

    /*
    @brief : Particles wiwhtin a cell are ordered accroding the flux values
    */
    void 
    ClemParticles::UpdateIndexAfterAdvection(){
        amrex::Print() << "Updating particle index after advection" << std::endl;
        const int lev = 0;

        for (MyPartIter pti(*this, lev); pti.isValid(); ++pti)
        {
            const int grid_id               = pti.index();
            const int tile_id               = pti.LocalTileIndex();

            auto& cell_fab                  = m_cell_vectors[grid_id];
            auto& particle_tile             = GetParticles(lev)[{grid_id, tile_id}];
            auto& particles                 = particle_tile.GetArrayOfStructs();
            const amrex::Box& bx            = pti.validbox();
            //amrex::Print() << "Processing grid: " << grid_id << " with box: " << bx << std::endl;
            for(amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)){
                auto& plist         = cell_fab(iv);
                int plist_size        = plist.size();
                
                auto mid_it         = std::stable_partition(plist.begin(), plist.end(),[&](int idx){ return particles[idx].idata(IntData::advected) == 0; });

                //Split the particle that already belong to the list adn the advected one.
                //the advected one separate by flux values first with the highest abd flux
                // if more than 1 cell are advected then order whith respect the index they have
                //Example p1.flux = 25 with idx 8 p2.flux = 25 with idx 7 p3.flux = 20 with idx 5
                //ORder     p1, p2, p3 since p1 and p2 have the same flux but p1 has higher index than p2

                //Update of the index of the particle according to the new order in the list
                //The standing particle may have an index p1 idx =3, p2 idx = 4 Now htey become 0 and 1 respectively then come the advected particles
                std::sort(plist.begin(), mid_it,[&](int a, int b) { return particles[a].idata(IntData::idx) < particles[b].idata(IntData::idx); });
                std::sort(mid_it, plist.end(), [&](int a, int b) {
                    const auto flux_a = particles[a].rdata(RealData::flux);
                    const auto flux_b = particles[b].rdata(RealData::flux);

                    if (flux_a == flux_b)
                        return particles[a].idata(IntData::idx) >
                               particles[b].idata(IntData::idx);

                    return flux_a > flux_b;
                });

                for (int n = 0; n < plist.size(); ++n){
                    particles[plist[n]].idata(IntData::idx) = n;
                    particles[plist[n]].idata(IntData::advected) = 0;
                    particles[plist[n]].rdata(RealData::flux) = 0.0;
                }
            }
        }
        m_cell_vector_updated = true;
    }

    /*
    @brief Initialized particles by "InitParticles" filled with the values from the Super cells
    */
    void
    ClemParticles::SetParticlePropertiesFromMF(const amrex::MultiFab& mf)
    {   
        amrex::MultiFab P(mf.boxArray(), mf.DistributionMap(), 1, 0 );
        clem::functionalities::ComputePressureFromMultiFab(P, mf);

        const int lev                   = 0;
        const auto& geom                = Geom(lev);
        const auto* dx                  = geom.CellSize();
        const amrex::Real vol           = (dx[0] * dx[0])/static_cast<amrex::Real>(NUM_LEM);

        AMREX_ASSERT(m_cell_vector_updated);

        for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const int grid_id       = mfi.index();
            const amrex::Box& bx    = mfi.validbox();
            const auto& mf_arr      = mf.array(mfi);
            const auto& P_arr       = P.array(mfi);
            auto& cell_fab          = m_cell_vectors[grid_id];

            // Loop over cells
            for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv))
            {
                amrex::Real massfrac[NUM_SPECIES] = {0.0};
                
                for(int sp = 0; sp < NUM_SPECIES ; sp++){
                    massfrac[sp] = mf_arr(iv, UFS + sp);
                }

                // Read cell data
                amrex::Real ueint       = mf_arr(iv, UEINT);
                amrex::Real rho         = mf_arr(iv, URHO);
                amrex::Real T           = mf_arr(iv, UTEMP);
                amrex::Real mass        = vol * rho;
                amrex::Real P           = P_arr(iv, 0);

                auto& p_indices         = cell_fab(iv);

                if (p_indices.empty()) continue;
                auto& particle_map = GetParticles(lev);
                for (auto& kv : particle_map)
                {
                    if (kv.first.first != grid_id) continue;

                    auto& particle_tile     = kv.second;
                    auto& particles         = particle_tile.GetArrayOfStructs();

                    for (int pidx : p_indices)
                    {   
                        //Filling PArticle with information
                        ParticleType& p             = particles[pidx];
                        p.rdata(RealData::mass)     = mass;
                        p.rdata(RealData::pressure) = P;
                        p.rdata(RealData::volume)   = vol;
                        p.rdata(RealData::flux)     = 0.0;
                        p.rdata(RealData::density)  = rho;
                        p.rdata(RealData::energy)   = ueint;
                        p.rdata(RealData::temp)     = T;

                        for(int sp = 0; sp < NUM_SPECIES ; sp++){
                            //amrex::Print() << RealData::SP(sp) << "  " << massfrac[sp] << std::endl;
                            p.rdata(RealData::SP(sp)) = massfrac[sp];
                        }
                    }
                }
            }
        }
        //WritePlotFile("plt_particles_values", "particles");
    }
    /*
    @brief The properties of the particles are comunicated to the super cell. Now need update since only captures the mass
    */
    void
    ClemParticles::GetCurrentAmountOfMassSpecies(amrex::MultiFab& mf)
    {
        const int lev = 0;
        const auto& geom                = Geom(lev);
        const auto* dx                  = geom.CellSize();
        const amrex::Real vol           = (dx[0] * dx[0]);

        AMREX_ASSERT(m_cell_vector_updated);
        mf.setVal(0.0);

        for (MyPartIter pti(*this, lev); pti.isValid(); ++pti)
        {
            const int grid_id = pti.index();
            const int tile_id = pti.LocalTileIndex();

            auto& particle_tile = GetParticles(lev)[{grid_id, tile_id}];
            auto& particles     = particle_tile.GetArrayOfStructs();

            // Loop over particles ONCE
            for (int pidx = 0; pidx < pti.numParticles(); ++pidx)
            {
                ParticleType& p = particles[pidx];
                amrex::IntVect iv(AMREX_D_DECL(p.idata(IntData::i),p.idata(IntData::j),p.idata(IntData::k)));

                amrex::Real mass_over_vol       = p.rdata(RealData::mass) / vol;
                amrex::Real one_rho             = 1.0 / p.rdata(RealData::density);

                for(int sp = 0; sp < NUM_SPECIES; sp++)
                    mf[pti](iv, sp) += p.rdata(RealData::SP(sp)) * mass_over_vol * one_rho;
            }
        }

        mf.SumBoundary(Geom(lev).periodicity());
    }

    void
    ClemParticles::ParticlesAdvection(amrex::Real dt)
    {
        //Here all celll has initially NUM_LEM particle.. this is enforced earlier.
        //PReviously  regriddign: all particle have the same dm and equal number of aprticles per cell
        const int lev                   = 0;
        const amrex::Geometry& geom     = Geom(lev);
        const amrex::Real* plo          = geom.ProbLo();
        const amrex::Real* dx           = geom.CellSize();
        const amrex::Box& domain        = Geom(lev).Domain();

        AMREX_ASSERT(m_cell_vector_updated);

        for(MyPartIter pti(*this, lev); pti.isValid(); ++pti){
            const int grid_id           = pti.index();
            const int tile_id           = pti.LocalTileIndex();

            auto& flux_arr              = FluxFaces[grid_id];
            auto& cell_fab              = m_cell_vectors[grid_id];

            const amrex::Box& bx        = pti.validbox();
            auto& particle_tile         = GetParticles(lev)[std::make_pair(grid_id, tile_id)];
            auto& particles             = particle_tile.GetArrayOfStructs();
            amrex::Vector<ParticleType> new_particles;

            for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv))
            {
                //amrex::Print() << "Processing cell: " << iv << std::endl;
                constexpr int numFaces                        = AMREX_SPACEDIM == 2 ? 4 : 6 ;
                advection::FluxDir fluxes[]         = { {flux_arr(iv, 0)* dt, advection::Face::Left},
                                                        {flux_arr(iv, 1)* dt, advection::Face::Right},
                                                        {flux_arr(iv, 2)* dt, advection::Face::Bottom},
                                                        {flux_arr(iv, 3)* dt, advection::Face::Top},
    #if AMREX_SPACEDIM == 3
                                                        {flux_arr(iv, 4)* dt, advection::Face::South},
                                                        {flux_arr(iv, 5)* dt, advection::Face::North},
    #endif                                            
                                                    };

                auto& p_indices             = cell_fab(iv);                                        
                int pointer_index           = NUM_LEM - 1;//We assume that all cell have NUM_LEM of particles

                //The fluxes are ordered
                //the highes abs flux to the smalles one. This is important since we want to move first the particle with the highest flux and then the one with the smaller flux. 
                //This is important since if we have more than 1 particle to move across a face we want to move first the one with the highest flux and then the one with the smaller flux.
                //Not that the flyu have sign. The positive flux will be advected to the next cell and the negative flux will be advected from the other side of the face.
                std::sort(std::begin(fluxes), std::end(fluxes),[](const advection::FluxDir& a, const advection::FluxDir& b) {return a.flux > b.flux;});

                for(int iface = 0; iface < numFaces; iface++){
                    if (pointer_index < 0) {
                        amrex::Print()<< "Minimin number of Lem achieved during splicing/Advection" << std::endl;
                        break;
                    }
                    //We only advect outwards the face with positive flux. The negative flux will be advected from the other side of the face
                    if (fluxes[iface].flux <= 0.0) continue;

                    const advection::Face face          = fluxes[iface].face;
                    const amrex::Real net_flux          = fluxes[iface].flux;
    

                    const amrex::IntVect dst     = advection::neighbor_iv(iv, face);
                    const amrex::IntVect src     = iv;

                    //To advect particle we only mode the particle to the destiantion. usually to the center of the next cell


                    amrex::Real flux_counter    = std::abs(net_flux);
                    //amrex::Print() << "Advection across face: " << static_cast<int>(face) << " with flux: " << net_flux << " index src: " << src << "  index des" << dst << std::endl;

                    while(flux_counter > 0.0 ){
                        //If the particle has not enough mass to complete the flux, then is fully advected
                        if( flux_counter > particles[p_indices[pointer_index]].rdata(RealData::mass)){
                            //amrex::Print() << "Moving entire particle" << std::endl;
                            ParticleType& p = particles[p_indices[pointer_index]];
                            //amrex::Print() << "Moving entire particle with mass: " << p.rdata(RealData::mass) << " Remaining flux to advect: " << flux_counter << std::endl;

                            //p.rdata(RealData::mass)     = flux_counter;
                            p.idata(IntData::advected)  = 1;
                            p.rdata(RealData::flux)     = std::abs(net_flux);

                            p.pos(0) =  plo[0] + (dst[0] + 0.5)*dx[0];
                            p.pos(1) =  plo[1] + (dst[1] + 0.5)*dx[1];
                            p.pos(2) =  plo[2] + (dst[2] + 0.5)*dx[2];
                            pointer_index  --;
                            flux_counter    -= p.rdata(RealData::mass);
                            //THe partile advected doe snot belong enymore to this cell so we remove it from the list of particle indices of this cell.
                            //After full advection the new destiantion cell will get it
                            p_indices.pop_back();
                           
                        }else{
                            //the mass is more thatn enough, Then we create a particle taht will be advected with the flux mass and the remaining mass will stay in the original particle
                            //amrex::Print() << "Splitting particles" << std::endl;
                            ParticleType p_out;
                            ParticleType& p_in      = particles[p_indices[pointer_index]];
                            p_out.id()              = ParticleType::NextID();
                            p_out.cpu()             = amrex::ParallelDescriptor::MyProc();

                            particle_functionalities::SetParticleWithOtherParticleValuesAdvection(iv, plo, dx, flux_counter, net_flux, p_in, p_out);
                            //const ParticleType p_adv  = CreateParticleAdvectionTypeDefault(dst, plo, dx, flux_counter, p_out);
                            p_in.rdata(RealData::mass)     -= flux_counter;
                            //amrex::Print() << "Mass of the original particle after splitting: " << p_in.rdata(RealData::mass) << " Mass of the new particle: " << p_out.rdata(RealData::mass) << std::endl;
                            new_particles.push_back(p_out);
                            flux_counter = 0.0;
                            break;
                        }
                    }
                }
            } //iv loop

            for (auto& p : new_particles) {
                particle_tile.push_back(p);
            }
        }//MPIter
        //flag to update the mapping grids and mesh
        m_cell_vector_updated = false;
    }


    void 
    ClemParticles::DoIsentropicPressureChange(amrex::MultiFab& Pressure)
    {   
        const int lev = 0;

        for (amrex::MFIter mfi(Pressure); mfi.isValid(); ++mfi)
        {
            const int grid_id       = mfi.index();
            const amrex::Box& bx    = mfi.validbox();
            const auto& P_arr       = Pressure.array(mfi);
            auto& cell_fab          = m_cell_vectors[grid_id];
            
            for(amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd();  bx.next(iv)){

                amrex::Real P               = P_arr(iv, 0);
                auto& plist                 = cell_fab(iv);
                int numberOfParticles       = plist.size();

                if (plist.empty()) continue;
                auto& particle_map = GetParticles(lev);

                for (auto& kv : particle_map)
                {
                    if (kv.first.first != grid_id) continue;
                    auto& particle_tile     = kv.second;
                    auto& particles         = particle_tile.GetArrayOfStructs();
                    for (int pidx : plist){  //Equivalent to iterate over all lem
                        ParticleType& p             = particles[pidx];
                        amrex::Real P_old           = p.rdata(RealData::pressure);

                        //Isentropic calculation
                        amrex::Real isentropicPChange   = std::pow(P/P_old, constants::isentropicRatio);
                        p.rdata(RealData::pressure)     = P;
                        p.rdata(RealData::temp)         *= isentropicPChange;
                    }
                }

            }
        }
    }//DoIsentropicPressureChange


    /*
    @brief This method creates a total number of particles euqla to NUM_LEM per cell.
    This is previous the DoRegridding function
    */
    void 
    ClemParticles::CreateParticleBeforeRegridding(){
        const int lev                   = 0;
        const auto& geom                = Geom(lev);
        const auto* dx                  = geom.CellSize();
        const auto* plo                 = geom.ProbLo();
        //     *  *  *          That why diviosn by N+1 lem element
        // -> --------- <-
        const amrex::Real dx_inner      = dx[0] / static_cast<amrex::Real>(NUM_LEM + 1);

        for(MyPartIter pti(*this, lev); pti.isValid(); ++pti){ 
            const int grid_id           = pti.index();
            const int tile_id           = pti.LocalTileIndex();
            auto& cell_fab              = m_cell_vectors[grid_id]; //Here the particles idex are stores
            const amrex::Box& bx        = pti.validbox();
            auto& particle_tile         = GetParticles(lev)[std::make_pair(grid_id, tile_id)];
            auto& particles             = particle_tile.GetArrayOfStructs();

            for(amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd();  bx.next(iv)){
                auto& plist         = cell_fab(iv);
                std::sort(plist.begin(), plist.end(),[&](int a, int b) {return particles[a].idata(IntData::idx)< particles[b].idata(IntData::idx);});
                
                const int numberOfLem   = plist.size();
                //Verification if the index are well sortex 
                for (int n = 0; n < numberOfLem; ++n) {
                    const int pidx = plist[n];
                    const int idx  = particles[pidx].idata(IntData::idx);
                    AMREX_ASSERT_WITH_MESSAGE(idx == n, "Particle indices are not sorted correctly before regridding");
                }

                if(numberOfLem < NUM_LEM){
                    amrex::Print() << "Creating new particles until NUM_LEM: " << NUM_LEM << std::endl;
                    //LEt say we have 5 particle and required 9. Then the last index to the aprticles is 4 and the last should be 8
                    for(int index_lem = numberOfLem ; index_lem < NUM_LEM; index_lem++){
                        ParticleType particle;
                        particle.id()      = ParticleType::NextID();
                        particle.cpu()     = amrex::ParallelDescriptor::MyProc();
                        particle_functionalities::SetParticleWithDefaultValues(iv, plo, dx, dx_inner, index_lem, particle);
                        //amrex::Print() << "Creating particle with index: " << index_lem << std::endl;
                        //ParticleType  p = CreateParticleTypeDefault(index_lem, iv, plo, dx, dx_inner);
                        particle_tile.push_back(particle);
                    }
                    //add particles until size = NUM_LEM
                }//we dont analyse when numberOfLem is equal or greter than NUM_LEM since will be deleted late

            }
        }
    }

    void
    ClemParticles::DoRegridding(){
        auto eos                    = pele::physics::PhysicsType::eos();

        CreateParticleBeforeRegridding();
        amrex::Print() << "Finish creating particles before regridding" << std::endl;
        Redistribute();
        amrex::Print() << "Finish Redistribute" << std::endl;
        DefineOwnershipParticleMesh();
        amrex::Print() << "Finish DefineOwnershipParticleMesh" << std::endl;

        //WritePlotFile("plt_Regridding_previous", "particles");
        const int lev                   = 0;
        const auto& geom                = Geom(lev);
        const auto* dx                  = geom.CellSize();
        const auto* plo                 = geom.ProbLo();
        //     *  *  *          That why diviosn by N+1 lem element
        // -> --------- <-
        const amrex::Real dx_inner      = dx[0] / static_cast<amrex::Real>(NUM_LEM + 1);

        AMREX_ASSERT(m_cell_vector_updated);

        for(MyPartIter pti(*this, lev); pti.isValid(); ++pti){ 
            const int grid_id           = pti.index();
            const int tile_id           = pti.LocalTileIndex();
            auto& cell_fab              = m_cell_vectors[grid_id]; //Here the particles idex are stores
            const amrex::Box& bx        = pti.validbox();
            auto& particle_tile         = GetParticles(lev)[std::make_pair(grid_id, tile_id)];
            auto& particles             = particle_tile.GetArrayOfStructs();

            //amrex::Print() << "Starting regridding procedure for cell: " << bx << std::endl;
            for(amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd();  bx.next(iv)){
                //amrex::Print() << "Regridding cell: " << iv << std::endl;
                auto& plist         = cell_fab(iv);
                std::sort(plist.begin(), plist.end(),[&](int a, int b) {return particles[a].idata(IntData::idx)< particles[b].idata(IntData::idx);});
                
                amrex::Real total_mass  = 0.0;
                const int numberOfLem   = plist.size();

                //Verification if the index are well sortex 
                for (int n = 0; n < numberOfLem; ++n) {
                    //if(n == 9) amrex::Print() << "Number of particle in the cell: " << numberOfLem << std::endl;
                    const int pidx = plist[n];
                    const int idx  = particles[pidx].idata(IntData::idx);
                    AMREX_ASSERT_WITH_MESSAGE(idx == n, "Particle indices are not sorted correctly before regridding");
                }


                for(const int pidx : plist){
                    total_mass += particles[pidx].rdata(RealData::mass);
                }
                const amrex::Real new_dm     = total_mass/static_cast<amrex::Real>(NUM_LEM);
                
                //Creation of a buffer with all the information of the particles before the regridding. This is required since we will be modifying the particles 
                //during the regridding and we need to keep the original information to do the mass average correctly
                //amrex::Print()<< iv << " Total mass: " << total_mass << " New dm: " << new_dm << "number of partice: " << numberOfLem << std::endl;

                constexpr int number_of_variable = NUM_SPECIES + 5; //5 is for mass, pressure, volume, density and temp
                amrex::Vector<amrex::Vector<amrex::Real, amrex::DefaultAllocator<amrex::Real>>> particles_data(numberOfLem, amrex::Vector<amrex::Real>(number_of_variable, 0.0));

                for(int i = 0; i < numberOfLem; i++){
                    ParticleType& p = particles[plist[i]];
                    particles_data[i][0] = p.rdata(RealData::mass);
                    particles_data[i][1] = p.rdata(RealData::pressure);
                    particles_data[i][2] = p.rdata(RealData::volume);
                    particles_data[i][3] = p.rdata(RealData::density);
                    particles_data[i][4] = p.rdata(RealData::temp);

                    amrex::Real one_rho = 1.0 / p.rdata(RealData::density);

                    for(int sp = 0; sp < NUM_SPECIES; sp++){
                        particles_data[i][5 + sp] = p.rdata(RealData::SP(sp)) * one_rho; 
                        //We store mass fraction since is more easy to do the mass average and then transform back to conservative after the regridding
                    }
                }

                //Set all particles values = 0.0 to avoid to use them during the regridding procedure. We will fill only until NUM_LEM, the rest of particles will be deleted at the end
                for(const int pidx : plist){
                    ParticleType& p = particles[pidx];
                    p.rdata(RealData::mass)     = 0.0;
                    p.rdata(RealData::pressure) = 0.0;
                    p.rdata(RealData::volume)   = 0.0;
                    p.rdata(RealData::density)  = 0.0;
                    p.rdata(RealData::temp)     = 0.0;
                    for(unsigned int sp = 0; sp < NUM_SPECIES ; sp++){
                        p.rdata(RealData::SP(sp)) = 0.0;
                    }
                }

                //Redistribution procedure

                int idx_source = 0;
                int idx_target = 0;
                while(idx_source < numberOfLem && idx_target < NUM_LEM){

                    ParticleType& p_target = particles[plist[idx_target]];

                    amrex::Real mass_source = particles_data[idx_source][0];
                    amrex::Real mass_need   = new_dm - p_target.rdata(RealData::mass);
                    //amrex::Print()<< "Source idx: " << idx_source << " Target idx: " << idx_target << " Mass source: " << mass_source << " Mass need: " << mass_need << std::endl;
                    if(mass_source >= mass_need && mass_need > 0.0){
                        //Store the sata in the particle target
                        //FIrst step just accumualte, alter we will divide by the new mass to have the mass average
                        //Mass fraction are stores in particles_data in primitive form
                        p_target.rdata(RealData::mass)     += mass_need;
                        p_target.rdata(RealData::pressure) += particles_data[idx_source][1]*mass_need;
                        p_target.rdata(RealData::volume)   += mass_need / particles_data[idx_source][3];
                        p_target.rdata(RealData::density)  = p_target.rdata(RealData::mass) / p_target.rdata(RealData::volume);
                        p_target.rdata(RealData::temp)     += particles_data[idx_source][4]*mass_need;

                        for(int sp = 0; sp < NUM_SPECIES; sp++){
                            p_target.rdata(RealData::SP(sp)) += particles_data[idx_source][5 + sp] * mass_need;
                        }

                        //Update the source particle data
                        particles_data[idx_source][0] -= mass_need; //mass
                        particles_data[idx_source][2] -= mass_need / particles_data[idx_source][3]; //volume
                        //Update the target particle data
                        idx_target++;
                    }else{
                        //Take all the mass from the source particle and update the target particle
                        p_target.rdata(RealData::mass)     += mass_source;
                        p_target.rdata(RealData::pressure) += particles_data[idx_source][1]*mass_source;
                        p_target.rdata(RealData::volume)   += mass_source / particles_data[idx_source][3];
                        p_target.rdata(RealData::density)  = p_target.rdata(RealData::mass) / p_target.rdata(RealData::volume);
                        p_target.rdata(RealData::temp)     += particles_data[idx_source][4]*mass_source;
                        for(int sp = 0; sp < NUM_SPECIES; sp++){
                            p_target.rdata(RealData::SP(sp)) += particles_data[idx_source][5 + sp]*mass_source;
                        }   
                        //Update the source particle data
                        particles_data[idx_source][0] = 0.0; //mass
                        particles_data[idx_source][2] = 0.0; //volume
                        //Update the target particle data
                        idx_source++;
                    }
                }

                //HEre we verify that the particle_data is empty
                for(int i = 0; i < numberOfLem; i++){
                    AMREX_ASSERT_WITH_MESSAGE(particles_data[i][0] == 0.0, "Not all the mass from the source particle has been used during the regridding procedure");
                }
                
                //Here we perfom mass average for the target particle, since during the regridding procedure we accumulate the mass 
                //and the other variables but we need to divide by the new mass to have the mass average. 
                //Note that if the mass is 0.0 then we dont do anything since the particle will be deleted later

                for(int i = 0; i < NUM_LEM; i++){
                    ParticleType& p = particles[plist[i]];

                    if(p.rdata(RealData::mass) > 0.0){
                        p.rdata(RealData::pressure) /= p.rdata(RealData::mass);
                        p.rdata(RealData::temp)     /= p.rdata(RealData::mass);
                        amrex::Real rho = 0.0;

                        for(int sp = 0; sp < NUM_SPECIES; sp++){
                            //The accumualtion is m_i*Y_i then to have the mass fraction we need to divide by the mass and then to
                            p.rdata(RealData::SP(sp)) /= p.rdata(RealData::volume);
                            rho += p.rdata(RealData::SP(sp));
                        }
                        p.rdata(RealData::density) = rho;
                    }
                }

                //needd to update ueint
                for(int i = 0; i < NUM_LEM; i++){
                    ParticleType& p                     = particles[plist[i]];
                    amrex::Real massfrac[NUM_SPECIES]   = {0.0};

                    amrex::Real rho         = p.rdata(RealData::density);
                    amrex::Real one_rho     = 1.0/p.rdata(RealData::density);
                    amrex::Real T           = p.rdata(RealData::temp);
                    amrex::Real energy      = 0.0;

                    for(int sp = 0; sp < NUM_SPECIES; sp++){
                        massfrac[sp] = p.rdata(RealData::SP(sp)) * one_rho;
                    }

                    eos.RTY2E(rho, T, massfrac, energy);
                    p.rdata(RealData::energy) = energy * p.rdata(RealData::density);
                }

                //eliminate particles with mass 0.0
                for(int i = 0; i < numberOfLem; i++){
                    ParticleType& p = particles[plist[i]];
                    //amrex::Print()<< "Particle " << i << " mass: " << p.rdata(RealData::mass) << std::endl;
                    if(p.rdata(RealData::mass) == 0.0){ 
                        //mark the particle for deletion by setting its index to a value higher than NUM_LEM, since we know that only particles with index lower than NUM_LEM are used during the regridding procedure
                        p.id() = -p.id(); //mark for deletion
                    }
                }
                //Reset the position and index of the stading aprticles
                for(int i = 0; i < NUM_LEM; i++){
                    ParticleType& p = particles[plist[i]];
                    p.pos(0)    = plo[0] + iv[0]*dx[0] + (i+1)*dx_inner;
                    p.pos(1)    = plo[1] + (iv[1] + 0.5)*dx[1];
#if (AMREX_SPACEDIM == 3)
                    p.pos(2)    = plo[2] + (iv[2] + 0.5)*dx[2];
#endif
                    p.idata(IntData::idx)      = i; 
                }

                for(int i = 0 ; i<NUM_LEM; i++){
                    ParticleType& p = particles[plist[i]];
                    //amrex::Print()<< "Particle " << i << " mass: " << p.rdata(RealData::mass) << std::endl;
                    if(p.rdata(RealData::mass) > 1.0){
                        amrex::Print()<< "Mass: " << p.rdata(RealData::mass) << " Pressure: " << p.rdata(RealData::pressure) << " Volume: " << p.rdata(RealData::volume) << " Density: " << p.rdata(RealData::density) << " Temp: " << p.rdata(RealData::temp) << std::endl;
                    }
                    AMREX_ASSERT_WITH_MESSAGE(p.rdata(RealData::mass) - new_dm < 1e-12, "Particle is not similar as the predicted new dm after the regridding procedure");
                    AMREX_ASSERT_WITH_MESSAGE(p.rdata(RealData::mass)  < 1, "mass cannot be more than 1");
                }
            }//bx iv
        }//mfi iter




        amrex::Print() << "Regridding procedure done for all the cells " << std::endl;
        Redistribute();
        amrex::Print() << "Resdistribute " << std::endl;
        //WritePlotFile("plt_Regridding_Intermediate", "particles");
        amrex::Print() << "Regridding procedure done for all the cells " << std::endl;
        

    }//DoRegridding

    void
    ClemParticles::DoClem(  pele::physics::transport::TransParm<pele::physics::PhysicsType::eos_type,pele::physics::PhysicsType::transport_type> const* ltransparm,
                            pele::physics::reactions::ReactorBase* reactor,
                            amrex::MultiFab& reactSrc, 
                            amrex::Real dt,
                         //Here variables to update superCell
                            amrex::MultiFab& nonReactSrc,
                            amrex::MultiFab& S_old,
                            amrex::MultiFab& S_new)
    {       
        const amrex::BoxArray& grids            = S_new.boxArray();
        const amrex::DistributionMapping& dmap  = S_new.DistributionMap();

        amrex::MultiFab PSuperGrid(grids, dmap, 1, 1);
        amrex::MultiFab NewAmountOfMassSpecies(grids, dmap, NUM_SPECIES, 1);
        amrex::MultiFab OldAmountOfMassSpecies(grids, dmap, NUM_SPECIES, 1);  
        amrex::MultiFab OmegaFiltered(grids, dmap,NUM_SPECIES, 1);

        NewAmountOfMassSpecies.setVal(0.0);
        OldAmountOfMassSpecies.setVal(0.0);
        OmegaFiltered.setVal(0.0);

        //Step 8 and 9
        //isentropic pressure change

        functionalities::ComputePressureFromMultiFab(PSuperGrid, S_new);
        DoIsentropicPressureChange(PSuperGrid);
        amrex::Print() << " ......... Isentropic pressure change" << std::endl;
        //WritePlotFile("plt_PIsentropic", "particles");

        amrex::ParallelDescriptor::Barrier();//-------------------------------------------------
        //step 10
        //Regridding
        Redistribute();
        amrex::Print() << " ......... Redistribute" << std::endl;
        DoRegridding();
        amrex::Print() << " ......... Do Regridding" << std::endl;
        Redistribute();
        amrex::Print() << " ......... Redistribute" << std::endl;
        DefineOwnershipParticleMesh();
        amrex::ParallelDescriptor::Barrier();//-------------------------------------------------
        amrex::Print() << " ......... Regridding" << std::endl;
        //WritePlotFile("plt_Regridding", "particles");
        
        //Step 11
        //Calculation or current omega
        GetCurrentAmountOfMassSpecies(OldAmountOfMassSpecies);
        amrex::Print() << " ......... Computation Old Mass Species" << std::endl;

        //step from 12 to 19
        const int lev                       = 0;
        const amrex::Geometry& geom         = Geom(lev);
        const amrex::Real* dx               = geom.CellSize();

        const amrex::Real areaFactor        = (AMREX_SPACEDIM == 2) ?  dx[0] : dx[0]*dx[0];
        /*
        Here Computation of:
            a) Mass diffusion
            b) Temperature Diffusion
            c) Reaction
        */
        for(MyPartIter pti(*this, lev); pti.isValid(); ++pti){
            const int grid_id           = pti.index();
            const int tile_id           = pti.LocalTileIndex();
            auto& cell_fab              = m_cell_vectors[grid_id]; //Here the particles idex are stores
            const amrex::Box& bx        = pti.validbox();
            auto& particle_tile         = GetParticles(lev)[std::make_pair(grid_id, tile_id)];
            auto& particles             = particle_tile.GetArrayOfStructs();
            
            for(amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd();  bx.next(iv)){
                auto& plist = cell_fab(iv);
                std::sort(plist.begin(), plist.end(),[&](int a, int b) {return particles[a].idata(IntData::idx)< particles[b].idata(IntData::idx);});
                
                //Create an array to store the data 
                std::array<amrex::Real, NUM_LEM> density{};
                std::array<amrex::Real, NUM_LEM> mass{};
                std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM> conservatives{}; //NUM_SPECIES + TEMP               
                std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM> primitives{}; //NUM_SPECIES + TEMP
                std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM> rhs{}; //NUM_SPECIES + TEMP
                std::array<std::array<amrex::Real, NUM_SPECIES + 1>, NUM_LEM> transport{}; //NUM_SPECIES + TEMP

                //Fill the the main data array with particle data.
                for(int lem = 0; lem < NUM_LEM; lem++){
                    ParticleType& p             = particles[plist[lem]];
                    particle_functionalities::CopyParticleDataToArray(mass, density, primitives, conservatives, rhs, p, lem);
                }

                //computes the diffusion and fill the transport array
                diffusion::ComputationOfTemperatureAndSpeciesDiffusion(mass, density, primitives, transport, rhs, areaFactor, ltransparm);

                amrex::Real current_time = 0.0;
                reactor->reactClemParticles(primitives, conservatives, rhs, density, dt, current_time
                #ifdef AMREX_USE_GPU
                        ,
                        amrex::Gpu::gpuStream()
                #endif
                ); //all conservtives container the new vall

                for(int lem = 0; lem < NUM_LEM; lem++){
                    ParticleType& particle             = particles[plist[lem]];
                    particle_functionalities::CopyArrayDataToParticle(mass, density, primitives, conservatives, rhs, particle, lem);            
                }
            }

        }
        //step 22
        GetCurrentAmountOfMassSpecies(NewAmountOfMassSpecies);
        amrex::Print() << " ......... Computation New Mass Species" << std::endl;
        functionalities::ComputeOmegaFiltered(NewAmountOfMassSpecies, OldAmountOfMassSpecies, OmegaFiltered, dt);
        amrex::Print() << " ......... Computation Omega Filtered" << std::endl;
        //Step 23
        ParticlesAdvection(dt);
        Redistribute();
        DefineOwnershipParticleMesh();
        //WritePlotFile("plt_AdvectionAfter", "particles");
        UpdateIndexAfterAdvection();
        //WritePlotFile("plt_Advection", "particles");
        amrex::Print() << " ......... Particle Advection" << std::endl;
        //Step 25
        functionalities::UpdatePeleCDataFromCLEM(OmegaFiltered, NewAmountOfMassSpecies, S_old, S_new, nonReactSrc, reactSrc, dt);
        amrex::Print() << " ......... Super-grid Update" << std::endl;
        amrex::Print() << "..... Reaction Finished" << std::endl;
        //WritePlotFile("plt_runing", "particles");
    }//DoClem,
}
