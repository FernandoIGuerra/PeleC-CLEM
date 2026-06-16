#ifndef CLEM_FUNCTIONALITIES_H
#define CLEM_FUNCTIONALITIES_H

#include <AMReX.H>
#include "PeleC.H"
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>

namespace clem {
    namespace functionalities{
        void ComputePressureFromMultiFab(amrex::MultiFab& P_supergrid, const amrex::MultiFab& mf);

        void ComputeOmegaFiltered( amrex::MultiFab& NewMassSpecies, 
                                                    amrex::MultiFab& OldMassSpecies, 
                                                    amrex::MultiFab& OmegaFiltered,
                                                    amrex::Real dt_les);
        
        void UpdatePeleCDataFromCLEM(amrex::MultiFab& omegaFiltered,
                                    amrex::MultiFab& newAmuntOfMassSpecies, //This the mass averge of all lem component
                                    amrex::MultiFab& S_old,
                                    amrex::MultiFab& S_new,
                                    amrex::MultiFab& nonReactSrc,
                                    amrex::MultiFab& reactSrc,
                                    amrex::Real dt_les);
    }
}

#endif /* CLEM_FUNCTIONALITIES_H */
