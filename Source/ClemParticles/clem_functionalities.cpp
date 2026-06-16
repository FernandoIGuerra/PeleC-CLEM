
#include "ClemParticles/clem_functionalities.H"
#include <array>

namespace clem {

    void functionalities::ComputePressureFromMultiFab(amrex::MultiFab& P_supergrid, const amrex::MultiFab& mf){
        for(amrex::MFIter mfi(mf); mfi.isValid(); ++mfi){
            const amrex::Box& bx                                = mfi.validbox();
            amrex::Array4<const amrex::Real> snew_array         = mf.array(mfi);
            const amrex::Array4<amrex::Real>& p_array           = P_supergrid.array(mfi);

            //Think if we have to pass Boxarray to react instead of the a unique lem element
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i,int j,int k) noexcept {
                auto eos = pele::physics::PhysicsType::eos();
                amrex::Real rho                     = snew_array(i,j,k, URHO);
                amrex::Real one_rho                 = 1.0/rho;
                amrex::Real T                       = snew_array(i,j,k, UTEMP);
                amrex::Real massloc[NUM_SPECIES]    = {0.0};
                amrex::Real eint                    = snew_array(i,j,k, UEINT) * one_rho;
                for(int sp = 0; sp < NUM_SPECIES; sp++){
                    massloc[sp] = snew_array(i,j,k, UFS + sp) * one_rho;
                } 
                amrex::Real P = 0.0;

                eos.RYET2P(rho, massloc, eint, T, P);
                //Here the Old and NewStates are update accordingly
                p_array(i,j,k) = P;

            });
        }
    }

    void functionalities::ComputeOmegaFiltered( amrex::MultiFab& NewMassSpecies, 
                                                amrex::MultiFab& OldMassSpecies, 
                                                amrex::MultiFab& OmegaFiltered,
                                                amrex::Real dt_les)
    {
        for(amrex::MFIter mfi(NewMassSpecies); mfi.isValid(); ++mfi){
        const amrex::Box& bx                        = mfi.validbox();
        const amrex::Array4<amrex::Real>& newY      = NewMassSpecies.array(mfi);
        const amrex::Array4<amrex::Real>& oldY      = OldMassSpecies.array(mfi);
        const amrex::Array4<amrex::Real>& omega     = OmegaFiltered.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept{
            for(int sp = 0; sp < NUM_SPECIES; sp++){
                omega(i,j,k, sp) = (newY(i,j,k, sp) - oldY(i,j,k, sp))/dt_les;
            }
            }); 
        }
    }

    void functionalities::UpdatePeleCDataFromCLEM   (amrex::MultiFab& omegaFiltered,
                                                    amrex::MultiFab& newAmuntOfMassSpecies, //This the mass averge of all lem component
                                                    amrex::MultiFab& S_old,
                                                    amrex::MultiFab& S_new,
                                                    amrex::MultiFab& nonReactSrc,
                                                    amrex::MultiFab& reactSrc,
                                                    amrex::Real dt_les)
    {
        for(amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi){
            const amrex::Box& bx                                = mfi.validbox();
            const amrex::Array4<amrex::Real>& omega             = omegaFiltered.array(mfi);
            const amrex::Array4<amrex::Real>& sold_arr          = S_old.array(mfi);
            const amrex::Array4<amrex::Real>& snew_arr          = S_new.array(mfi);
            const amrex::Array4<amrex::Real>& nonrs_arr         = nonReactSrc.array(mfi);
            const amrex::Array4<amrex::Real>& newY_array        = newAmuntOfMassSpecies.array(mfi);
            auto const& I_R                                     = reactSrc.array(mfi);

            // CLEM coupling: LES (supergrid) owns rho, rho*u, rho*E_hydro. The
            // continuity equation is solved exclusively on the LES grid, so density
            // and momentum here MUST NOT be overwritten -- otherwise the acoustic
            // mode carried by hydro is destroyed.
            //
            // The LEM/subgrid contributes only:
            //   (a) composition (mass fractions), projected onto hydro density,
            //   (b) chemistry heat release into the LES energy equation, where
            //       qdot = -sum_i h_i * omega[sp]  and  omega is the pure LEM
            //       (diffusion + reaction) rate already computed in DoClem.
            //
            // I_R is exposed as the LEM subgrid source rate (not mixed with the
            // LES non-reactive source nonrs_arr).
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept{
                const amrex::Real rho_hydro = snew_arr(i, j, k, URHO);

                // Project LEM composition onto hydro density. Sum to rho_hydro by
                // construction, so mass-fraction normalization is exact.
                amrex::Real rho_lem_sum = 0.0;
                for (int sp = 0; sp < NUM_SPECIES; sp++)
                    rho_lem_sum += newY_array(i, j, k, sp);

                const amrex::Real lem_inv = (rho_lem_sum > 0.0) ? 1.0 / rho_lem_sum : 0.0;
                for (int sp = 0; sp < NUM_SPECIES; sp++)
                    snew_arr(i, j, k, UFS + sp) = rho_hydro * newY_array(i, j, k, sp) * lem_inv;

                // LEM-only source rate for diagnostics / AMR reflux.
                for (int sp = 0; sp < NUM_SPECIES; sp++)
                    I_R(i, j, k, sp) = omega(i, j, k, sp);
            });

            // Heat release from LEM contribution. Add to UEINT and UEDEN; do not
            // recompute kinetic energy (momentum is unchanged from hydro).
            //Here we include the Sub grid effect on the whole model
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                auto eos = pele::physics::PhysicsType::eos();

                amrex::Real hi[NUM_SPECIES]    = {0.0};
                amrex::Real Yspec[NUM_SPECIES] = {0.0};

                const amrex::Real rho_hydro = snew_arr(i, j, k, URHO);
                const amrex::Real rho_inv   = 1.0 / rho_hydro;
                for (int sp = 0; sp < NUM_SPECIES; sp++)
                    Yspec[sp] = snew_arr(i, j, k, UFS + sp) * rho_inv;

                eos.RTY2Hi(rho_hydro, snew_arr(i, j, k, UTEMP), Yspec, hi);

                amrex::Real qdot = 0.0;
                for (int sp = 0; sp < NUM_SPECIES; sp++)
                    qdot -= hi[sp] * omega(i, j, k, sp);

                snew_arr(i, j, k, UEINT) += dt_les * qdot;
                snew_arr(i, j, k, UEDEN) += dt_les * qdot;

                I_R(i, j, k, NUM_SPECIES)     = qdot;
                I_R(i, j, k, NUM_SPECIES + 1) = qdot;
            });
        }
    }
}