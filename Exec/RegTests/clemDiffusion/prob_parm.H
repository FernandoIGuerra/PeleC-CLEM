#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real p_init = 1013250.0; // 1 atm
  amrex::Real Y_H2_C = 0.06;
  amrex::Real Y_O2_C = 0.5;
  amrex::Real Y_N2_C = 0.44;
  amrex::Real T_C = 940.0;
  amrex::Real rho_C = 0.0;
  amrex::Real e_C = 0.0;

  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac_C = {0.0};
  amrex::Real Y_H2_H = 0.06;
  amrex::Real Y_O2_H = 0.5;
  amrex::Real Y_N2_H = 0.44;
  amrex::Real T_H = 940.0;
  amrex::Real rho_H = 0.0;
  amrex::Real e_H = 0.0;

  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac_H = {0.0};

  amrex::Real xliml = 0.0;
  amrex::Real xlimh = 0.0;
};

struct ProbParmHost
{
  ProbParmHost() = default;
};
#endif
