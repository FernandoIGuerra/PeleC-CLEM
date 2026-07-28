#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

// ============================================================
// 1D freely-propagating premixed flame for CLEM accuracy tests.
//
// The whole domain is initialized by interpolating a precomputed
// 1D flame (PMF / Tecplot POINT format) onto the mesh, so the flame
// starts at its correct thickness, structure and velocity field and
// only fast acoustic transients need to relax.  The cold-boundary
// state of that profile (with u = S_L) is also used as the Hard
// inflow at x-lo.  Outflow (FOExtrap) at x-hi.
//
// Reference: phi = 0.4 H2/air, LiDryer, p = 1 atm, Tu = 298 K
//            (LiDryer_H2_p1_phi0_4000tu0300.dat, S_L ~ 22.8 cm/s)
//
// All units are CGS (PeleC convention).
// ============================================================
struct ProbParmDevice
{
  amrex::Real pamb = 1013250.0; // 1 atm in dyn/cm^2

  // Flame placement: the reference flame is queried at (x - standoff).
  // standoff is auto-set in amrex_probinit so the reference flame front
  // lands at x = x_flame in this domain; it can be overridden in inputs.
  amrex::Real x_flame = 0.35; // desired front position in this domain [cm]
  amrex::Real standoff = 0.0; // filled in amrex_probinit (or queried)

  // --- PMF reference profile (filled by read_pmf) ---
  int pmf_N = 0;             // number of points
  int pmf_M = 0;             // number of variables (temp, u, rho, species...)
  bool pmf_do_average = false;
  amrex::Real* d_pmf_X = nullptr;
  amrex::Real* d_pmf_Y = nullptr;
};

struct ProbParmHost
{
  amrex::Vector<amrex::Real> h_pmf_X;
  amrex::Vector<amrex::Real> h_pmf_Y;
  amrex::Gpu::DeviceVector<amrex::Real> pmf_X;
  amrex::Gpu::DeviceVector<amrex::Real> pmf_Y;
  ProbParmHost() : pmf_X(0), pmf_Y(0) {}
};

#endif
