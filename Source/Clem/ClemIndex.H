#ifndef CLEMINDEX_H
#define CLEMINDEX_H

#include <AMReX_REAL.H>
#include "mechanism.H"

// Fernando-Clem: compile-time number of LEM elements per LES cell
// (alpaca-style: override at build time with -DCLEM_NLEM=<n>)
#ifndef CLEM_NLEM
#define CLEM_NLEM 12
#endif

// Fernando-Clem: compile-time splicing policy for the strategy comparison
//   0 = side_aware (Maxwell steps 23-24: donor from the exit-facing end,
//       arrivals spliced onto the end matching their entry face)
//   1 = append_downstream (legacy: donor from the back, all arrivals appended
//       at the downstream end - known artificial front diffusivity)
#ifndef CLEM_SPLICE_POLICY
#define CLEM_SPLICE_POLICY 0
#endif

namespace clem {

namespace config {
// Fernando-Clem: compile-time configuration, accessed as constexpr everywhere
static constexpr int n_lem = CLEM_NLEM;

enum class SplicePolicy : int { side_aware = 0, append_downstream = 1 };
static constexpr SplicePolicy splice_policy =
  static_cast<SplicePolicy>(CLEM_SPLICE_POLICY);
} // namespace config

// Fernando-Clem: per-particle real components p_k = [rho, m, v, T, eint, p,
// Y_k] (primitive storage, Maxwell's element state). press is the element's
// own pressure: reference of the isentropic LES->SGS stage. xflux is
// transient splicing bookkeeping: the face mass flux that carried the
// parcel, used to order same-end arrivals; cleared on reindex
struct RealData
{
  static constexpr int rho = 0;   // density
  static constexpr int mass = 1;  // element mass
  static constexpr int vol = 2;   // element volume
  static constexpr int T = 3;     // temperature
  static constexpr int eint = 4;  // specific internal energy
  static constexpr int press = 5; // element pressure
  static constexpr int xflux = 6; // |face mass flux| of the last splice move
  static constexpr int Y0 = 7;    // first species mass fraction, Y_k = Y0 + k
  static constexpr int ncomps = 7 + NUM_SPECIES;
};

// Fernando-Clem: per-particle int components. entry_side is transient
// splicing bookkeeping (0 = resident, 1 = entered through the left end,
// 2 = entered right, 3 = transverse arrival); cleared on reindex
struct IntData
{
  static constexpr int lem_index = 0;  // position along the 1-D LEM line
  static constexpr int entry_side = 1; // splicing entry code, see above
  static constexpr int ncomps = 2;
};

} // namespace clem
#endif
