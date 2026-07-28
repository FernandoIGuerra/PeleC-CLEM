#ifndef CLEMPARTICLECONTAINER_H
#define CLEMPARTICLECONTAINER_H

#include <AMReX_AmrParticles.H>

#include "ClemIndex.H"

namespace clem {

// Fernando-Clem: iterators mirror the spray pattern (MyParIter/MyParConstIter)
class ClemParIter
  : public amrex::ParIter<RealData::ncomps, IntData::ncomps, 0, 0>
{
public:
  using amrex::ParIter<RealData::ncomps, IntData::ncomps, 0, 0>::ParIter;
};

class ClemParConstIter
  : public amrex::ParConstIter<RealData::ncomps, IntData::ncomps, 0, 0>
{
public:
  using amrex::ParConstIter<RealData::ncomps, IntData::ncomps, 0, 0>::
    ParConstIter;
};

// Fernando-Clem: storage-only container (AoS, no SoA comps, like spray);
// physics operators live outside so strategies can be swapped at compile time
class ClemParticleContainer
  : public amrex::AmrParticleContainer<RealData::ncomps, IntData::ncomps, 0, 0>
{
public:
  using BaseType =
    amrex::AmrParticleContainer<RealData::ncomps, IntData::ncomps, 0, 0>;
  using ParticleType = BaseType::ParticleType;

  explicit ClemParticleContainer(amrex::AmrCore* amr) : BaseType(amr) {}

  // Fernando-Clem: single-level ctor so the container is testable without a
  // full Amr hierarchy
  ClemParticleContainer(
    const amrex::Geometry& geom,
    const amrex::DistributionMapping& dmap,
    const amrex::BoxArray& ba)
    : BaseType({geom}, {dmap}, {ba}, {})
  {
  }
};

} // namespace clem
#endif
