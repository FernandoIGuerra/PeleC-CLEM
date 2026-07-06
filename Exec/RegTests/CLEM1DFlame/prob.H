#ifndef PROB_H
#define PROB_H

#include <AMReX_ParmParse.H>
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

// Linear interpolation of the reference flame profile at coordinate [xlo, xhi].
// y_vector[0] = temperature, [1] = velocity, [2] = density,
// y_vector[3 + n] = species mole fraction.  (Mirrors Exec/RegTests/PMF.)
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
pmf(
  const amrex::Real xlo,
  const amrex::Real xhi,
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4>& y_vector,
  const ProbParmDevice& prob_parm)
{
  amrex::Real xmid = 0.5 * (xlo + xhi);
  int loside = -1;
  int hiside = -1;
  if (xmid < prob_parm.d_pmf_X[0]) {
    loside = 0;
    hiside = 0;
  }
  if (xmid > prob_parm.d_pmf_X[prob_parm.pmf_N - 1]) {
    loside = prob_parm.pmf_N - 1;
    hiside = prob_parm.pmf_N - 1;
  }
  if (loside == -1) {
    for (int i = 0; i < prob_parm.pmf_N - 1; i++) {
      if (
        (xmid >= prob_parm.d_pmf_X[i]) && (xmid <= prob_parm.d_pmf_X[i + 1])) {
        loside = i;
        hiside = i + 1;
      }
    }
  }
  for (int j = 0; j < prob_parm.pmf_M; j++) {
    const amrex::Real x1 = prob_parm.d_pmf_X[loside];
    const amrex::Real y1 = prob_parm.d_pmf_Y[prob_parm.pmf_N * j + loside];
    const amrex::Real x2 = prob_parm.d_pmf_X[hiside];
    const amrex::Real y2 = prob_parm.d_pmf_Y[prob_parm.pmf_N * j + hiside];
    amrex::Real dydx = (loside == hiside) ? 0.0 : (y2 - y1) / (x2 - x1);
    y_vector[j] = y1 + dydx * (xmid - x1);
  }
}

// Fill the conserved state from the reference profile at file-coordinate
// (x - standoff), with the flame propagating in x (velocity into UMX).
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
fill_state_from_pmf(
  const amrex::Real xq, amrex::Real* s, ProbParmDevice const& prob_parm)
{
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {{0.0}};
  pmf(xq, xq, pmf_vals, prob_parm);

  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  for (int n = 0; n < NUM_SPECIES; n++) {
    molefrac[n] = pmf_vals[3 + n];
  }
  const amrex::Real T = pmf_vals[0];
  const amrex::Real u = pmf_vals[1]; // streamwise velocity (x)

  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  auto eos = pele::physics::PhysicsType::eos();
  eos.X2Y(molefrac, massfrac);
  amrex::Real rho, energy;
  eos.PYT2RE(prob_parm.pamb, massfrac, T, rho, energy);

  s[URHO] = rho;
  s[UMX] = rho * u;
  s[UMY] = 0.0;
  s[UMZ] = 0.0;
  s[UEINT] = rho * energy;
  s[UEDEN] = rho * (energy + 0.5 * u * u);
  s[UTEMP] = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    s[UFS + n] = rho * massfrac[n];
  }
}

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

  amrex::Real s[NVAR] = {0.0};
  fill_state_from_pmf(x - prob_parm.standoff, s, prob_parm);
  for (int n = 0; n < NVAR; n++) {
    state(i, j, k, n) = s[n];
  }
}

// Inflow (x-lo, "Hard"): cold premixed reactants from the cold end of the
// reference profile (its velocity there is S_L).  x-hi is outflow (FOExtrap).
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const amrex::Real* /*s_int[NVAR]*/,
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real /*time*/,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*turb_fluc*/)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  // sgn = +1 -> low side, sgn = -1 -> high side
  const amrex::Real xb = (sgn == 1) ? prob_lo[idir] : prob_hi[idir];
  fill_state_from_pmf(xb - prob_parm.standoff, s_ext, prob_parm);
}

void pc_prob_close();

using ProblemSpecificFunctions = DefaultProblemSpecificFunctions;

#endif
