#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

// using namespace amrex::literals;

struct ProbParmDevice
{

  //Channel
  amrex::Real Pres_amb = 1013250.0;
  amrex::Real Temp_amb = 1000.0;
  amrex::Real Yox_amb = 0.23;

  amrex::Real Yox_channel = 0.23;
  amrex::Real Mach_channel = 1.5;

  // amrex::Real vel_channel = 930.0*100;

  //JET

  amrex::Real Pres_jet = 1.29 * 1013250.0; // 1.29 atm in Ba
  amrex::Real temp_jet = 300.0;
  amrex::Real Yfuel_jet = 1.0;
  amrex::Real Mach_jet = 1.02;

  amrex::Real vel_jet = 132100.0;
  amrex::Real massflow_jet = 0.109; //g/s
  amrex::Real global_equivalence_ratio = 0.15;
  amrex::Real jet_to_freestream_momentum_ratio = 0.58;

  // C1 case
  amrex::Real centx = 4.2;

  amrex::Real centz = 0.25;
  amrex::Real r_hole = 0.05; // d = 1mm


  int inject_fuel = 0;  //MODIFICAR




  //int init_type = 0;
  //amrex::Real cavity_depth = 0.305;
};

struct ProbParmHost
{
  ProbParmHost() = default;
};

#endif
