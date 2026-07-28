#ifndef CLEMSCRATCH_H
#define CLEMSCRATCH_H

#include <AMReX_Vector.H>

#include "ClemIndex.H"

namespace clem {

// Fernando-Clem: per-rank reusable scratch arena for the LEM line operators.
//
// A line operator cannot work on the particle AoS in place: it needs the cell's
// elements gathered in LINE order (the map hands them out in AoS order), it
// needs derived state that is not stored on the particle (X_k, h_k, rho*D_k,
// rho*lambda, Cp), and it needs face and rate buffers that have no particle
// counterpart at all. That is a lot of copied information per cell, and sizing
// it per cell would put n_cell * n_substep heap allocations inside the inner
// loop of every step.
//
// So the buffers are owned by ONE arena per rank, grown on the first cell of
// the run and reused unchanged for every cell of every step afterwards (resize
// never shrinks). The arena is thread_local, so if a line operator is ever put
// inside an OpenMP region each thread gets its own line and the buffers stay
// race-free by construction; with MPI-only there is exactly one arena per rank.
//
// Layout: element quantities are indexed [l], per-species quantities are
// species-major per element, [l * NUM_SPECIES + k]; face quantities are indexed
// by the interior face f, which sits between elements f and f+1 (so face f is
// the "l+1/2" face of element l = f). The LEM line is PERIODIC (Maxwell's
// homogeneous sub-grid sample): face n_elem-1 wraps element n_elem-1 -> 0, so
// there are n_elem faces. Inter-cell transport stays the splicing operator's job.
struct LineScratch
{
  // Fernando-Clem: the per-rank instance; all line operators share it
  static LineScratch& get()
  {
    static thread_local LineScratch arena;
    return arena;
  }

  // Fernando-Clem: gathered element state (line order, lem_index ascending)
  amrex::Vector<int> pidx;          // index of the element in the tile AoS
  amrex::Vector<amrex::Real> mass;  // element mass = the mass-coordinate dm
  amrex::Vector<amrex::Real> rho;
  amrex::Vector<amrex::Real> T;
  amrex::Vector<amrex::Real> eint;
  amrex::Vector<amrex::Real> press; // element pressure, held constant over the
                                    // sub-cycle (the LEM array is isobaric)
  amrex::Vector<amrex::Real> Y;

  // Fernando-Clem: derived state, recomputed from (rho, T, Y) every sub-step so
  // the driving forces and the coefficients stay consistent as the line evolves
  amrex::Vector<amrex::Real> X;      // mole fractions
  amrex::Vector<amrex::Real> hk;     // species enthalpies (formation included)
  amrex::Vector<amrex::Real> rhoD;   // rho * (rho D_k W_k/W) - the mass-coordinate
                                     // species coefficient
  amrex::Vector<amrex::Real> rholam; // rho * lambda - the mass-coordinate
                                     // thermal coefficient
  amrex::Vector<amrex::Real> cp;     // mixture Cp
  amrex::Vector<amrex::Real> cpk;    // per-species Cp_k; the temperature
                                     // equation's interspecies term uses Cp_k,
                                     // NOT h_k (that would be the energy form)
  amrex::Vector<amrex::Real> wbar;   // mean molecular weight; the species
                                     // stability limit needs W/W_k (see
                                     // DiffusionOperator::stableSubStep)

  // Fernando-Clem: interior-face fluxes
  amrex::Vector<amrex::Real> Sk;   // corrected species mass flux (mass equation)
  amrex::Vector<amrex::Real> qf;   // conductive heat flux (rho lambda) dT/dm
  amrex::Vector<amrex::Real> qint; // per-face interspecies product
                                   // sum_k Cp_k S_k dT/dm (temperature equation)

  // Fernando-Clem: rates of the current sub-step
  amrex::Vector<amrex::Real> dYdt;
  amrex::Vector<amrex::Real> dTdt; // T is the transported thermal variable

  int n_elem = 0; // elements on the gathered line
  int n_face = 0; // faces = n_elem (periodic line: face n_elem-1 wraps to 0)

  // Fernando-Clem: grow to hold a line of n elements. Steady state is reached
  // on the first cell, so every later call is just the two assignments
  void resize(const int n)
  {
    n_elem = n;
    n_face = n; // periodic: one face per element (the last wraps n-1 -> 0)

    if (static_cast<int>(mass.size()) >= n) {
      return; // Fernando-Clem: already big enough - the common case
    }

    const int ns = n * NUM_SPECIES;
    const int nf = n; // periodic line

    pidx.resize(n);
    mass.resize(n);
    rho.resize(n);
    T.resize(n);
    eint.resize(n);
    press.resize(n);
    cp.resize(n);
    wbar.resize(n);
    rholam.resize(n);
    dTdt.resize(n);

    Y.resize(ns);
    X.resize(ns);
    hk.resize(ns);
    cpk.resize(ns);
    rhoD.resize(ns);
    dYdt.resize(ns);

    qf.resize(nf);
    qint.resize(nf);
    Sk.resize(nf * NUM_SPECIES);
  }
};

} // namespace clem
#endif
