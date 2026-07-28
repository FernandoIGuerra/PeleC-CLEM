#ifndef CLEMREGRID_H
#define CLEMREGRID_H

#include "ClemManager.H"

namespace clem {

// Fernando-Clem: regridding operator (Maxwell thesis, p. 172). Rebins each
// LES cell's LEM line into exactly config::n_lem equal-mass elements while
// conserving mass, volume, species mass and internal energy. Stateless
// operator class (works through the manager's public API) so alternative
// regridding strategies can be swapped at compile time.
class Regridder
{
public:
  static void regrid(ClemManager& mgr, int lev);

private:
  // Fernando-Clem: top up cells owning fewer than n_lem elements with
  // zero-mass placeholders; they become merge targets and inherit mass
  static void ensureCellCapacity(ClemManager& mgr, int lev);

  // Fernando-Clem: per-cell two-pointer source->target equal-mass merge.
  // Note the particle state is primitive (T, specific e, mass fractions Y_k),
  // so the merge is a direct mass-weighted average - none of the old
  // conservative (rho Y, rho e) conversions are needed
  static void mergeCells(ClemManager& mgr, int lev);
};

} // namespace clem
#endif
