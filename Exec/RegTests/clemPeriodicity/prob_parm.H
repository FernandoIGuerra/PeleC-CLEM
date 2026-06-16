#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real p_init = 1013250.0; // 1 atm
  amrex::Real Y_H2_low = 0.06;
  amrex::Real Y_O2_low = 0.5;
  amrex::Real Y_N2_low = 0.44;
  amrex::Real T_low = 940.0;
  amrex::Real rho_low = 0.0;
  amrex::Real e_low = 0.0;

  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac_low = {0.0};
  amrex::Real Y_H2_high = 0.06;
  amrex::Real Y_O2_high = 0.5;
  amrex::Real Y_N2_high = 0.44;
  amrex::Real T_high = 940.0;
  amrex::Real rho_high = 0.0;
  amrex::Real e_high = 0.0;

  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac_high = {0.0};

  amrex::Real yliml = 0.0;
  amrex::Real ylimh = 0.0;

  amrex::Real velx_low = 5.0;
  amrex::Real vely_low = 0.0;

  amrex::Real velx_high = -5.0;
  amrex::Real vely_high = 0.0;
};

struct ProbParmHost
{
  ProbParmHost() = default;
};
#endif
