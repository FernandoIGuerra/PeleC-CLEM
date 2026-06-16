#ifndef PROB_H
#define PROB_H

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

#include "mechanism.H"

#include "PeleC.H"
#include "IndexDefines.H"
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
  // Set the state
  if(x>prob_parm.xliml && x<prob_parm.xlimh){ //Region to diffuse
    state(i, j, k, URHO) = prob_parm.rho_H;
    state(i, j, k, UMX) = 0.0;
    state(i, j, k, UMY) = 0.0;
    state(i, j, k, UMZ) = 0.0;
    state(i, j, k, UEINT) = prob_parm.rho_H * prob_parm.e_H;
    state(i, j, k, UEDEN) = prob_parm.rho_H * prob_parm.e_H;
    state(i, j, k, UTEMP) = prob_parm.T_H;
    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i, j, k, UFS + n) = prob_parm.rho_H * prob_parm.massfrac_H[n];
    }
  }else{
    state(i, j, k, URHO) = prob_parm.rho_C;
    state(i, j, k, UMX) = 0.0;
    state(i, j, k, UMY) = 0.0;
    state(i, j, k, UMZ) = 0.0;
    state(i, j, k, UEINT) = prob_parm.rho_C * prob_parm.e_C;
    state(i, j, k, UEDEN) = prob_parm.rho_C * prob_parm.e_C;
    state(i, j, k, UTEMP) = prob_parm.T_C;
    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i, j, k, UFS + n) = prob_parm.rho_C * prob_parm.massfrac_C[n];
    }
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const amrex::Real* /*s_int[NVAR]*/,
  amrex::Real* /*s_ext[NVAR]*/,
  const int /*idir*/,
  const int /*sgn*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& /*prob_parm*/,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*turb_fluc*/)
{
}

void pc_prob_close();

using ProblemSpecificFunctions = DefaultProblemSpecificFunctions;

#endif
