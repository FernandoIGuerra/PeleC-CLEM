#ifndef PROB_H
#define PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

#include "mechanism.H"

#include "PeleC.H"
#include "IndexDefines.H"
#include "Constants.H"
#include "PelePhysics.H"
#include "Tagging.H"
#include "ProblemSpecificFunctions.H"
#include "prob_parm.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_initdata(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  //Here we implementa a shear flowin the y direction
  // Set the state

  if( y<prob_parm.ylimh){ //Region to diffuse
    state(i, j, k, URHO) = prob_parm.rho_low;
    state(i, j, k, UMX) = prob_parm.velx_low * prob_parm.rho_low;
    state(i, j, k, UMY) = prob_parm.vely_low * prob_parm.rho_low;
    state(i, j, k, UMZ) = 0.0;
    state(i, j, k, UEINT) = prob_parm.rho_low * ( prob_parm.e_low + 0.5 * 
                                                ( prob_parm.velx_low * prob_parm.velx_low + 
                                                  prob_parm.vely_low * prob_parm.vely_low));

    state(i, j, k, UEDEN) = prob_parm.rho_low * prob_parm.e_low;
    state(i, j, k, UTEMP) = prob_parm.T_low;
    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i, j, k, UFS + n) = prob_parm.rho_low * prob_parm.massfrac_low[n];
    }
  }else{
    state(i, j, k, URHO) = prob_parm.rho_high;
    state(i, j, k, UMX) = prob_parm.velx_high * prob_parm.rho_high;
    state(i, j, k, UMY) = prob_parm.vely_high * prob_parm.rho_high;
    state(i, j, k, UMZ) = 0.0;
    state(i, j, k, UEINT) = prob_parm.rho_high * (prob_parm.e_high + 0.5 * 
                                                (prob_parm.velx_high * prob_parm.velx_high + 
                                                 prob_parm.vely_high * prob_parm.vely_high));

    state(i, j, k, UEDEN) = prob_parm.rho_high * prob_parm.e_high;
    state(i, j, k, UTEMP) = prob_parm.T_high;
    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i, j, k, UFS + n) = prob_parm.rho_high * prob_parm.massfrac_high[n];
    }
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const amrex::Real* /*s_int[NVAR]*/,
  amrex::Real s_ext[NVAR],
  const int /*idir*/,
  const int /*sgn*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*turb_fluc*/)
{


    s_ext[URHO] = prob_parm.rho_low;
    s_ext[UMX] = prob_parm.velx_low * prob_parm.rho_low;
    s_ext[UMY] = prob_parm.vely_low * prob_parm.rho_low;
    s_ext[UMZ] = 0.0;
    s_ext[UEDEN] = prob_parm.rho_low * ( prob_parm.e_low + 0.5 * 
                                                ( prob_parm.velx_low * prob_parm.velx_low + 
                                                  prob_parm.vely_low * prob_parm.vely_low));
                                                  
    s_ext[UEINT] = prob_parm.rho_low * prob_parm.e_low;
    s_ext[UTEMP] = prob_parm.T_low;
    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[UFS + n] = prob_parm.rho_low * prob_parm.massfrac_low[n];
    }

}

void pc_prob_close();

using ProblemSpecificFunctions = DefaultProblemSpecificFunctions;

#endif
