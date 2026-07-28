#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

// ============================================================
// CLEMTest2D: 2D jet-in-crossflow reactor.
//
// A cold H2 jet is injected normal (+y) into a hot air crossflow (+x)
// through a slot in the y-lo wall. Downstream of the slot the jet plume
// bends over, mixes with the crossflow across a shear layer, and (given
// enough residence time/temperature) ignites - the classic "jet in
// crossflow" (JICF) configuration, used here as a genuine 2D CLEM reactor
// from which spatial mixing/ignition statistics can be extracted.
//
// This case deliberately puts CLEM's subgrid line (1-D along x only, no
// transverse structure) in a flow with strong transverse (y) bulk motion:
// the jet's own splicing hits the same "transverse arrival" degradation
// documented for solid-body rotation (Maxwell limitation #2). CLEMTest2D
// exists to quantify that degradation in a real 2D reacting flow, not to
// dodge it - see memory clem-transport-validation.
//
// All units are CGS (PeleC convention).
// ============================================================
struct ProbParmDevice
{
  amrex::Real pamb = 1013250.0; // ambient pressure, 1 atm [dyn/cm^2]

  // --- Crossflow (x-lo Hard inflow; also the t=0 fill of the whole domain)
  amrex::Real T_cross = 1200.0;  // K - hot enough for H2/air to autoignite
  amrex::Real u_cross = 1500.0;  // cm/s (+x), ~15 m/s duct velocity
  amrex::Real Y_O2_cross = 0.232; // standard air mass fractions
  amrex::Real Y_N2_cross = 0.768;

  // --- Jet (slot in the y-lo wall; rest of that wall is a no-slip wall)
  amrex::Real T_jet = 300.0;    // K, cold pure H2
  amrex::Real v_jet = 6000.0;   // cm/s (+y), ~60 m/s -> jet/crossflow
                                 // velocity ratio ~4 (typical JICF range 2-10)
  amrex::Real jet_xc = 0.3;     // slot center [cm]
  amrex::Real jet_width = 0.1;  // slot width [cm]
  amrex::Real jet_xlo = 0.0;    // filled in amrex_probinit
  amrex::Real jet_xhi = 0.0;    // filled in amrex_probinit
};

struct ProbParmHost
{
  ProbParmHost() = default;
};

#endif
