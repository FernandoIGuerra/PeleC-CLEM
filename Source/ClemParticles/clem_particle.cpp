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

    void
    ClemParticles::WritePlotFileNamed(const std::string& dir, const std::string& name) const
    {
        amrex::Vector<std::string> real_names = {
            "mass", "pressure", "volume", "flux", "density", "energy", "temp"
        };
        amrex::Vector<std::string> spec_names;
        pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(spec_names);
        for (const auto& sn : spec_names) {
            real_names.push_back("Y_" + sn);
        }
        const amrex::Vector<std::string> int_names = {
            "idx", "i", "j", "k", "advected"
        };
        WritePlotFile(dir, name, real_names, int_names);
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
    ClemParticles::DefineOwnershipParticleMesh(bool update_particles_after_advection)
    {

        const int lev = 0;
        //We will evaluate later when should we avoid this
        //if (m_cell_vector_updated) return;

        //We clean all the ownership relations
        for(amrex::MFIter mfi = MakeMFIter(lev, false); mfi.isValid(); ++mfi){
            const int grid_id           = mfi.index();
            auto& fab = m_cell_vectors[grid_id];
            fab.resize(mfi.validbox());
            fab.setVal(std::vector<int>{});
        }

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
                amrex::IntVect iv           = this->Index(p, lev); //Get the index from the position i = (x/dx) 
                m_cell_vectors[grid_id](iv).push_back(pindex); //add to each cell the belonging cells
            }

        }
        m_cell_vector_updated = true;
        //This is a general algorthmin that sort the particles according to their index and flux 

        if (update_particles_after_advection){   
            UpdateIndexAfterAdvection();}
    }

    /*
    @brief : Particles wiwhtin a cell are ordered accroding the flux values
    */
    void
    ClemParticles::UpdateIndexAfterAdvection(){

        const int lev                   = 0;
        const auto& geom                = Geom(lev);
        const auto* dx                  = geom.CellSize();
        const auto* plo                 = geom.ProbLo();
        const amrex::Real dx_inner      = dx[0] / static_cast<amrex::Real>(NUM_LEM + 1);

        for (MyPartIter pti(*this, lev); pti.isValid(); ++pti)
        {
            const int grid_id               = pti.index();
            const int tile_id               = pti.LocalTileIndex();
            auto& cell_fab                  = m_cell_vectors[grid_id];
            auto& particle_tile             = GetParticles(lev)[{grid_id, tile_id}];
            auto& particles                 = particle_tile.GetArrayOfStructs();
            const amrex::Box& bx            = pti.validbox();


            for(amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)){
                auto& plist         = cell_fab(iv);
                int plist_size        = plist.size();
                
                auto mid_it         = std::stable_partition(plist.begin(), plist.end(),[&](int idx){ return particles[idx].idata(IntData::advected) == 0; });

                //Split the particle that already belong to the list and the advected one.
                //the advected one separate by flux values first with the highest flux
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
                    ParticleType& p                  = particles[plist[n]];
                    p.idata(IntData::idx)            = n;
                    p.idata(IntData::advected)       = 0;
                    p.rdata(RealData::flux)          = 0.0;

                    //Reset the particle position within the cell according to its new LEM index,
                    //following the same convention as DoRegridding. Only do this for the NUM_LEM
                    //particles that will remain after the upcoming regridding - repositioning the
                    //extra ones (n >= NUM_LEM, to be merged/deleted) could push pos() into a
                    //neighboring cell and desync it from idata(i/j/k)/cell_fab.
                    if(n < NUM_LEM){
                        AMREX_D_TERM(
                            p.pos(0) = plo[0] + iv[0]*dx[0] + (n+1)*dx_inner;,
                            p.pos(1) = plo[1] + (iv[1] + 0.5)*dx[1];,
                            p.pos(2) = plo[2] + (iv[2] + 0.5)*dx[2];
                        )
                    }
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
            const int tile_id       = mfi.LocalTileIndex();
            const amrex::Box& bx    = mfi.validbox();
            const auto& mf_arr      = mf.array(mfi);
            const auto& P_arr       = P.array(mfi);
            auto& cell_fab          = m_cell_vectors[grid_id];
            auto& particle_tile     = GetParticles(lev)[std::make_pair(grid_id, tile_id)];
            auto& particles         = particle_tile.GetArrayOfStructs();

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
                        p.rdata(RealData::SP(sp)) = massfrac[sp];
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
    ClemParticles::ParticlesAdvection(amrex::Real dt, const amrex::MultiFab& state_bc)
    {
        //Here all celll has initially NUM_LEM particle.. this is enforced earlier.
        //PReviously  regriddign: all particle have the same dm and equal number of aprticles per cell
        const int lev                   = 0;
        const amrex::Geometry& geom     = Geom(lev);
        const amrex::Real* plo          = geom.ProbLo();
        const amrex::Real* phi          = geom.ProbHi();
        const amrex::Real* dx           = geom.CellSize();
        const amrex::Box& domain        = Geom(lev).Domain();
        const auto& periodicity         = geom.periodicity();
        //inflow implementation - EOS used to recover the pressure of the incoming fluid from its (rho, T, Y) ghost-cell state
        auto eos                        = pele::physics::PhysicsType::eos();
        
        // Create periodic boundary condition flags
        amrex::IntVect is_periodic{AMREX_D_DECL(
                                                    periodicity.isPeriodic(0) ? 1 : 0,
                                                    periodicity.isPeriodic(1) ? 1 : 0,
                                                    periodicity.isPeriodic(2) ? 1 : 0
                                                )};

        //outflow implementation - physical BC types per direction (PCPhysBCType encoding, outflow == 2), set from PeleC::phys_bc via SetPhysBC
        const amrex::IntVect& lo_bc = m_lo_bc;
        const amrex::IntVect& hi_bc = m_hi_bc;

        for(MyPartIter pti(*this, lev); pti.isValid(); ++pti){
            const int grid_id           = pti.index();
            const int tile_id           = pti.LocalTileIndex();
            auto& flux_arr              = FluxFaces[grid_id];
            auto& cell_fab              = m_cell_vectors[grid_id];
            const auto state_arr        = state_bc[grid_id].const_array();  //inflow implementation - ghost cells of this grid hold the bcnormal inflow Dirichlet state
            const amrex::Box& bx        = pti.validbox();
            auto& particle_tile         = GetParticles(lev)[std::make_pair(grid_id, tile_id)];
            auto& particles             = particle_tile.GetArrayOfStructs();

            amrex::Vector<ParticleType> new_particles;

            for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv))
            {   
                // CALCULATION OF THE FLUXES
                constexpr int numFaces                          = AMREX_SPACEDIM == 2 ? 4 : 6 ;
                advection::FluxDir fluxes[]                     = { {flux_arr(iv, 0)* dt, advection::Face::Left},
                                                                    {flux_arr(iv, 1)* dt, advection::Face::Right},
                                                                    {flux_arr(iv, 2)* dt, advection::Face::Bottom},
                                                                    {flux_arr(iv, 3)* dt, advection::Face::Top},
                #if AMREX_SPACEDIM == 3
                                                                    {flux_arr(iv, 4)* dt, advection::Face::South},
                                                                    {flux_arr(iv, 5)* dt, advection::Face::North},
                #endif                                            
                                                                };

                auto& p_indices             = cell_fab(iv);                                        
                int pointer_index           = NUM_LEM - 1;
                //This must be true since previous this step regridding is done
                //We assume that all cell have NUM_LEM of particles

                //The fluxes are ordered
                //the highes abs flux to the smalles one. This is important since we want to move first the particle with the highest flux and then the one with the smaller flux. 
                //This is important since if we have more than 1 particle to move across a face we want to move first the one with the highest flux and then the one with the smaller flux.
                //Not that the flyu have sign. The positive flux will be advected to the next cell and the negative flux will be advected from the other side of the face.
                std::sort(std::begin(fluxes), std::end(fluxes),[](const advection::FluxDir& a, const advection::FluxDir& b) {return a.flux > b.flux;});

                for(int iface = 0; iface < numFaces; iface++){

                    const advection::Face face          = fluxes[iface].face;
                    const amrex::Real net_flux          = fluxes[iface].flux;
                    const amrex::IntVect dst            = advection::neighbor_iv(iv, face);
                    const amrex::IntVect src            = iv;

                    //We only advect outwards across faces with negative (outgoing) flux. 
                    //Positive (incoming) flux is supplied by the neighbour cell as ITS outflow.

                    //INFLOW BOUNDARY CONDITIONS: For the inflow check is the clel is a next to the BC
                    if (net_flux >= 0.0) {
                        //inflow implementation - ...except at a domain inflow boundary, where there is no interior neighbour to provide the mass.
                        //There the mass flows from OUTSIDE in: the source is the ghost cell (outside the domain) and the destination is the
                        //boundary cell iv (inside). We create a new particle in iv carrying the incoming mass with the LES ghost-cell (bcnormal) fluid state.
                        //Guard: dst outside domain ⟹ iv is the adjacent boundary cell.
                        //Interior cells always have dst inside the domain so they can never inject inflow particles.
                        if (net_flux > 0.0 && !domain.contains(dst) && boundaryCondition::is_inflow_boundary_face(dst, domain, is_periodic, lo_bc, hi_bc))
                        {
                            //the outward neighbour across an inflow boundary face is the ghost cell that holds the incoming-fluid state
                            const amrex::IntVect ghost_src = dst;

                            ParticleType p_inflow;
                            p_inflow.id()  = ParticleType::NextID();
                            p_inflow.cpu() = amrex::ParallelDescriptor::MyProc();

                            //incoming-fluid thermodynamic state is read from the ghost cell filled by the LES boundary condition (bcnormal)
                            const amrex::Real rho_in   = state_arr(ghost_src, URHO);
                            const amrex::Real one_rho  = 1.0 / rho_in;
                            const amrex::Real T_in     = state_arr(ghost_src, UTEMP);

                            amrex::Real massfrac[NUM_SPECIES] = {0.0};
                            for (int sp = 0; sp < NUM_SPECIES; ++sp) {
                                massfrac[sp] = state_arr(ghost_src, UFS + sp) * one_rho;
                            }
                            amrex::Real P_in = 0.0;
                            eos.RTY2P(rho_in, T_in, massfrac, P_in);

                            //the particle enters the boundary cell iv (inside the domain); place it at the cell centre so ownership maps to iv
                            AMREX_D_TERM(
                                p_inflow.pos(0) = plo[0] + (iv[0] + 0.5) * dx[0];,
                                p_inflow.pos(1) = plo[1] + (iv[1] + 0.5) * dx[1];,
                                p_inflow.pos(2) = plo[2] + (iv[2] + 0.5) * dx[2];
                            )
                            AMREX_D_TERM(
                                p_inflow.idata(IntData::i) = iv[0];,
                                p_inflow.idata(IntData::j) = iv[1];,
                                p_inflow.idata(IntData::k) = iv[2];
                            )
                            p_inflow.idata(IntData::idx)      = 100; //reassigned by UpdateIndexAfterAdvection
                            p_inflow.idata(IntData::advected) = 1;

                            p_inflow.rdata(RealData::mass)     = net_flux;          //mass entering through the face over dt
                            p_inflow.rdata(RealData::density)  = rho_in;
                            p_inflow.rdata(RealData::volume)   = net_flux * one_rho;
                            p_inflow.rdata(RealData::temp)     = T_in;
                            p_inflow.rdata(RealData::pressure) = P_in;
                            p_inflow.rdata(RealData::energy)   = state_arr(ghost_src, UEINT);  //rho*eint, conserved
                            p_inflow.rdata(RealData::flux)     = net_flux;
                            for (int sp = 0; sp < NUM_SPECIES; ++sp) {
                                p_inflow.rdata(RealData::SP(sp)) = state_arr(ghost_src, UFS + sp);  //rho*Y, conserved
                            }
                            new_particles.push_back(p_inflow);
                        }
                        continue;
                    }

                    if (pointer_index < 0) {
                        amrex::Print()<< "Minimin number of Lem achieved during splicing/Advection" << std::endl;
                        break;
                    }

                    //outflow implementation - mass crossing this face leaves the subgrid domain through a non-periodic outflow boundary
                    const bool exits_through_outflow = boundaryCondition::is_outside_domain_on_outflow(dst, domain, is_periodic, lo_bc, hi_bc);

                    //To advect particle we only mode the particle to the destiantion. usually to the center of the next cell


                    amrex::Real flux_counter    = std::abs(net_flux);

                    while(flux_counter > 0.0 ){
                        //If the particle has not enough mass to complete the flux, then is fully advected
                        if( flux_counter > particles[p_indices[pointer_index]].rdata(RealData::mass)){
                            //amrex::Print() << "Moving entire particle" << std::endl;
                            ParticleType& p = particles[p_indices[pointer_index]];
                            //amrex::Print() << "Moving entire particle with mass: " << p.rdata(RealData::mass) << " Remaining flux to advect: " << flux_counter << std::endl;

                            if (exits_through_outflow) {
                                //outflow implementation - the whole particle leaves the domain; mark it for deletion so Redistribute() purges it
                                flux_counter -= p.rdata(RealData::mass);
                                p.id() = -1;
                            } else {
                                //p.rdata(RealData::mass)     = flux_counter;
                                p.idata(IntData::advected)  = 1;
                                p.rdata(RealData::flux)     = std::abs(net_flux);

                                // Apply periodic BC wrapping to destination cell
                                auto [dst_wrapped, reentered] = boundaryCondition::wrap_iv(dst, domain, is_periodic);
                                //since dst_wrapped is already in the correct form of the doamin then no further calculation are required to set the new position of the particle
                                p.pos(0) =  plo[0] + (dst_wrapped[0] + 0.5)*dx[0];
                                p.pos(1) =  plo[1] + (dst_wrapped[1] + 0.5)*dx[1];
                                p.pos(2) =  plo[2] + (dst_wrapped[2] + 0.5)*dx[2];

                                //Update the cell indices to match the new cell, otherwise GetCurrentAmountOfMassSpecies
                                //would keep scattering this particle's mass/species into its old cell
                                AMREX_D_TERM(
                                    p.idata(IntData::i) = dst_wrapped[0];,
                                    p.idata(IntData::j) = dst_wrapped[1];,
                                    p.idata(IntData::k) = dst_wrapped[2];
                                )

                                flux_counter    -= p.rdata(RealData::mass);
                            }

                            pointer_index  --;
                            //THe partile advected doe snot belong enymore to this cell so we remove it from the list of particle indices of this cell.
                            //After full advection the new destiantion cell will get it
                            p_indices.pop_back();

                        }else{
                            ParticleType& p_in      = particles[p_indices[pointer_index]];

                            if (exits_through_outflow) {
                                //outflow implementation - only the split-off mass leaves the domain; no new particle is created, p_in stays with the remainder
                                p_in.rdata(RealData::mass)     -= flux_counter;
                                flux_counter = 0.0;
                                break;
                            }

                            //the mass is more thatn enough, Then we create a particle taht will be advected with the flux mass and the remaining mass will stay in the original particle
                            //amrex::Print() << "Splitting particles" << std::endl;
                            ParticleType p_out;
                            p_out.id()              = ParticleType::NextID();
                            p_out.cpu()             = amrex::ParallelDescriptor::MyProc();

                            // For split particles, set destination cell with periodic BC wrapping
                            auto [dst_wrapped, reentered] = boundaryCondition::wrap_iv(dst, domain, is_periodic);
                            // Set particle properties from source particle
                            //here take the dst_wrapped to set the positions
                            particle_functionalities::SetParticleWithOtherParticleValuesAdvection(dst_wrapped, plo, dx, flux_counter, net_flux, p_in, p_out);



                            p_in.rdata(RealData::mass)     -= flux_counter;
                            new_particles.push_back(p_out);
                            flux_counter = 0.0;
                            break;
                        }
                    }
                } //iterate over faces
            } //iv loop

            //ADD ALL NEW PARTICLES
            for (auto& p : new_particles) {
                particle_tile.push_back(p);
            }

        }//MPIter
    }

    void
    ClemParticles::DoIsentropicPressureChange(amrex::MultiFab& Pressure)
    {   
        const int lev = 0;

        for (amrex::MFIter mfi(Pressure); mfi.isValid(); ++mfi)
        {
            const int grid_id       = mfi.index();
            const int tile_id       = mfi.LocalTileIndex();
            const amrex::Box& bx    = mfi.validbox();
            const auto& P_arr       = Pressure.array(mfi);
            auto& cell_fab          = m_cell_vectors[grid_id];
            auto& particle_tile     = GetParticles(lev)[std::make_pair(grid_id, tile_id)];
            auto& particles         = particle_tile.GetArrayOfStructs();

            for(amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd();  bx.next(iv)){

                amrex::Real P               = P_arr(iv, 0);
                auto& plist                 = cell_fab(iv);

                if (plist.empty()) continue;

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
    }//DoIsentropicPressureChange


    /**
     * @brief This method creates a total number of particles euqla to NUM_LEM per cell. This is previous the DoRegridding function
     */
    void 
    ClemParticles::CreateParticleBeforeRegridding(){
        const int lev                   = 0;
        const auto& geom                = Geom(lev);
        const auto* dx                  = geom.CellSize();
        const auto* plo                 = geom.ProbLo();
        // locate the particle
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

        //Fill all LES cells with NUM_LEM particles
        CreateParticleBeforeRegridding();
        //Particles are distributed between ranks
        Redistribute();
        //We defined 
        DefineOwnershipParticleMesh();

        //WritePlotFile("plt_Regridding_previous", "particles");
        const int lev                   = 0;
        const auto& geom                = Geom(lev);
        const auto* dx                  = geom.CellSize();
        const auto* plo                 = geom.ProbLo();

        const amrex::Real dx_inner      = dx[0] / static_cast<amrex::Real>(NUM_LEM + 1);

        for(MyPartIter pti(*this, lev); pti.isValid(); ++pti){ 
            const int grid_id           = pti.index();
            const int tile_id           = pti.LocalTileIndex();
            auto& cell_fab              = m_cell_vectors[grid_id]; //Here the particles idex are stores
            const amrex::Box& bx        = pti.validbox();
            auto& particle_tile         = GetParticles(lev)[std::make_pair(grid_id, tile_id)];
            auto& particles             = particle_tile.GetArrayOfStructs();

            constexpr int number_of_variable = NUM_SPECIES + 5; //5 is for mass, pressure, volume, density and temp
            amrex::Vector<amrex::Real> particles_data;          //Scratch buffer reused across cells via resize() to avoid per-cell heap allocations


            for(amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd();  bx.next(iv)){
                auto& plist         = cell_fab(iv);
                std::sort(plist.begin(), plist.end(),[&](int a, int b) {return particles[a].idata(IntData::idx)< particles[b].idata(IntData::idx);});
                
                amrex::Real total_mass  = 0.0;
                const int numberOfLem   = plist.size();

                //Verification if the index are well sortex 
                for (int n = 0; n < numberOfLem; ++n) {
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

                particles_data.resize(numberOfLem * number_of_variable);

                for(int i = 0; i < numberOfLem; i++){
                    ParticleType& p = particles[plist[i]];
                    particles_data[i*number_of_variable + 0] = p.rdata(RealData::mass);
                    particles_data[i*number_of_variable + 1] = p.rdata(RealData::pressure);
                    particles_data[i*number_of_variable + 2] = p.rdata(RealData::volume);
                    particles_data[i*number_of_variable + 3] = p.rdata(RealData::density);
                    particles_data[i*number_of_variable + 4] = p.rdata(RealData::temp);

                    //avoid 1/0 = inf for placeholder particles (density == 0), which would turn their
                    //zero species data into NaN and poison the redistribution below even when mass_source == 0
                    amrex::Real one_rho = (p.rdata(RealData::density) > 0.0) ? 1.0 / p.rdata(RealData::density) : 0.0;

                    for(int sp = 0; sp < NUM_SPECIES; sp++){
                        particles_data[i*number_of_variable + 5 + sp] = p.rdata(RealData::SP(sp)) * one_rho;
                        //We store mass fraction since is more easy to do the mass average and then transform back to conservative after the regridding
                    }
                }

                //Set all particles values = 0.0 to avoid to use them during the regridding procedure. 
                //We will fill only until NUM_LEM, the rest of particles will be deleted at the end
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

                    amrex::Real mass_source         = particles_data[idx_source*number_of_variable + 0];
                    amrex::Real mass_need           = new_dm - p_target.rdata(RealData::mass);
                    const amrex::Real src_density   = particles_data[idx_source*number_of_variable + 3];

                    if(mass_source >= mass_need && mass_need > 0.0){
                        //Store the sata in the particle target
                        //First step just accumualte, alter we will divide by the new mass to have the mass average
                        //Mass fraction are stores in particles_data in primitive form
                        const amrex::Real vol_need = (src_density > 0.0) ? mass_need / src_density : 0.0; //avoid 0/0 from zero-density placeholder particles
                        p_target.rdata(RealData::mass)     += mass_need;
                        p_target.rdata(RealData::pressure) += particles_data[idx_source*number_of_variable + 1] * mass_need;
                        p_target.rdata(RealData::volume)   += vol_need;
                        p_target.rdata(RealData::density)  = (p_target.rdata(RealData::volume) > 0.0) ? p_target.rdata(RealData::mass) / p_target.rdata(RealData::volume) : 0.0;
                        p_target.rdata(RealData::temp)     += particles_data[idx_source*number_of_variable + 4] * mass_need;

                        for(int sp = 0; sp < NUM_SPECIES; sp++){
                            p_target.rdata(RealData::SP(sp)) += particles_data[idx_source*number_of_variable + 5 + sp] * mass_need;
                        }

                        //Update the source particle data
                        particles_data[idx_source*number_of_variable + 0] -= mass_need; //mass
                        particles_data[idx_source*number_of_variable + 2] -= vol_need; //volume
                        //Update the target particle data
                        idx_target++;
                    }else{
                        //Take all the mass from the source particle and update the target particle
                        const amrex::Real vol_source = (src_density > 0.0) ? mass_source / src_density : 0.0; //avoid 0/0 from zero-density placeholder particles
                        p_target.rdata(RealData::mass)     += mass_source;
                        p_target.rdata(RealData::pressure) += particles_data[idx_source*number_of_variable + 1]*mass_source;
                        p_target.rdata(RealData::volume)   += vol_source;
                        p_target.rdata(RealData::density)  = (p_target.rdata(RealData::volume) > 0.0) ? p_target.rdata(RealData::mass) / p_target.rdata(RealData::volume) : 0.0;
                        p_target.rdata(RealData::temp)     += particles_data[idx_source*number_of_variable + 4]*mass_source;
                        for(int sp = 0; sp < NUM_SPECIES; sp++){
                            p_target.rdata(RealData::SP(sp)) += particles_data[idx_source*number_of_variable + 5 + sp]*mass_source;
                        }
                        //Update the source particle data
                        particles_data[idx_source*number_of_variable + 0] = 0.0; //mass
                        particles_data[idx_source*number_of_variable + 2] = 0.0; //volume
                        //Update the target particle data
                        idx_source++;
                    }
                }

                //Here we verify that the particle_data is empty
                constexpr amrex::Real eps = 1e-15; 
                for (int i = 0; i < numberOfLem; ++i) {
                    amrex::Real val = particles_data[i * number_of_variable + 0];
                    AMREX_ASSERT_WITH_MESSAGE(std::abs(val) < eps,"Mass not fully consumed; residual = " << val);
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
                    if(p.rdata(RealData::mass) == 0.0) continue; // placeholder; will be deleted below
                    amrex::Real massfrac[NUM_SPECIES]   = {0.0};

                    amrex::Real rho         = p.rdata(RealData::density);
                    amrex::Real one_rho     = (rho > 0.0) ? 1.0/rho : 0.0;
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
                    AMREX_ASSERT_WITH_MESSAGE(p.rdata(RealData::mass) - new_dm < 1e-12, "Particle is not similar as the predicted new dm after the regridding procedure");
                }
            }//bx iv
        }//mfi iter
    }//DoRegridding

    void
    ClemParticles::DoClem(  pele::physics::transport::TransParm<pele::physics::PhysicsType::eos_type,pele::physics::PhysicsType::transport_type> const* ltransparm,
                            pele::physics::reactions::ReactorBase* reactor,
                            amrex::MultiFab& reactSrc, 
                            amrex::Real dt,
                         //Here variables to update superCell
                            amrex::MultiFab& nonReactSrc,
                            amrex::MultiFab& S_old,
                            amrex::MultiFab& S_new,
                         //inflow implementation - FillPatched state with bcnormal-filled ghost cells for incoming-fluid properties
                            const amrex::MultiFab& Sborder)
    {
        /**
         * S_new = U^** <= 0.5(u^n + u^*) + 0.5(S^{n+1}(u^*) + I_R) dt
         * S_old = u^n
         * */    
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

        amrex::ParallelDescriptor::Barrier();
        //step 10
        //Regridding
        Redistribute();
        DoRegridding();
        Redistribute();
        DefineOwnershipParticleMesh();
        WritePlotFileNamed("plt_AfterRegridding", "particles");
        amrex::ParallelDescriptor::Barrier();
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
                diffusion::ComputationOfTemperatureAndSpeciesDiffusion(mass, density, primitives, transport, rhs, areaFactor, ltransparm, dt);

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

                //Diagnostic: dump the post-reaction LEM state for the suspected cell so the
                //evolution of mass/species/temp can be tracked across the steps where the
                //instability develops. Remove once the root cause is identified.
                if (iv == amrex::IntVect(AMREX_D_DECL(0,29,0))) {
                    for(int lem = 0; lem < NUM_LEM; lem++){
                        ParticleType& p = particles[plist[lem]];
                        amrex::Print() << "[CLEM-debug] cell " << iv
                                        << " lem " << lem
                                        << " idx="      << p.idata(IntData::idx)
                                        << " mass="     << p.rdata(RealData::mass)
                                        << " density="  << p.rdata(RealData::density)
                                        << " volume="   << p.rdata(RealData::volume)
                                        << " temp="     << p.rdata(RealData::temp)
                                        << " pressure=" << p.rdata(RealData::pressure);
                        for(int sp = 0; sp < NUM_SPECIES; sp++){
                            amrex::Print() << " Y" << sp << "=" << p.rdata(RealData::SP(sp));
                        }
                        amrex::Print() << std::endl;
                    }
                }
            }

        }
        //step 22
        GetCurrentAmountOfMassSpecies(NewAmountOfMassSpecies);
        functionalities::ComputeOmegaFiltered(NewAmountOfMassSpecies, OldAmountOfMassSpecies, OmegaFiltered, dt);

        //Step 23
        //outflow implementation - ParticlesAdvection marks particles that leave through outflow boundaries (id = -1); Redistribute() then purges them
        //inflow implementation - ParticlesAdvection injects new particles at inflow boundaries using Sborder's ghost-cell (bcnormal) state
        ParticlesAdvection(dt, Sborder);
        Redistribute();
        DefineOwnershipParticleMesh();
        UpdateIndexAfterAdvection();

        //Step 25
        functionalities::UpdatePeleCDataFromCLEM(OmegaFiltered, NewAmountOfMassSpecies, S_old, S_new, nonReactSrc, reactSrc, dt);
        WritePlotFileNamed("plt_AfterDoClem", "particles");
    }//DoClem,
}
