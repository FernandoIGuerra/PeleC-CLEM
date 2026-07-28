#ifndef CLEMADVECTION_H
#define CLEMADVECTION_H

#include <utility>

#include <AMReX.H>
#include <AMReX_IntVect.H>
#include <AMReX_Box.H>

// Fernando-Clem: face/flux bookkeeping and boundary helpers for the splicing
// operator, ported from the old clem_advection.H

namespace clem {

namespace advection {

enum class Face { Left = 0, Right = 1, Bottom = 2, Top = 3, South = 4, North = 5 };

static constexpr int num_faces = 2 * AMREX_SPACEDIM;

struct FluxDir
{
  amrex::Real flux;
  Face face;
};

constexpr amrex::IntVect face_offsets[] = {
  amrex::IntVect{AMREX_D_DECL(-1, 0, 0)}, // Left
  amrex::IntVect{AMREX_D_DECL(1, 0, 0)},  // Right
  amrex::IntVect{AMREX_D_DECL(0, -1, 0)}, // Bottom
  amrex::IntVect{AMREX_D_DECL(0, 1, 0)},  // Top
#if AMREX_SPACEDIM == 3
  amrex::IntVect{AMREX_D_DECL(0, 0, -1)}, // South
  amrex::IntVect{AMREX_D_DECL(0, 0, 1)},  // North
#endif
};

AMREX_FORCE_INLINE
amrex::IntVect
neighbor_iv(const amrex::IntVect& iv, Face face)
{
  return iv + face_offsets[static_cast<int>(face)];
}

} // namespace advection

namespace boundary {

// Fernando-Clem: wrap a destination cell through periodic boundaries;
// second member reports whether wrapping occurred
AMREX_FORCE_INLINE
std::pair<amrex::IntVect, bool>
wrap_iv(
  const amrex::IntVect& iv,
  const amrex::Box& domain,
  const amrex::IntVect& is_periodic)
{
  amrex::IntVect r = iv;
  if (domain.contains(r)) {
    return {r, false};
  }
  bool was_wrapped = false;
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    if (is_periodic[dir] != 0) {
      const int length = domain.length(dir);
      if (r[dir] < domain.smallEnd(dir)) {
        r[dir] += length;
        was_wrapped = true;
      } else if (r[dir] > domain.bigEnd(dir)) {
        r[dir] -= length;
        was_wrapped = true;
      }
    }
  }
  return {r, was_wrapped};
}

// Fernando-Clem: cell outside the domain across a non-periodic outflow
// boundary (PCPhysBCType encoding: outflow == 2)
AMREX_FORCE_INLINE
bool
is_outside_domain_on_outflow(
  const amrex::IntVect& iv,
  const amrex::Box& domain,
  const amrex::IntVect& is_periodic,
  const amrex::IntVect& lo_bc,
  const amrex::IntVect& hi_bc)
{
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    if (is_periodic[dir] != 0) {
      continue;
    }
    if (iv[dir] < domain.smallEnd(dir) && lo_bc[dir] == 2) {
      return true;
    }
    if (iv[dir] > domain.bigEnd(dir) && hi_bc[dir] == 2) {
      return true;
    }
  }
  return false;
}

// Fernando-Clem: crossing this face enters the domain through a non-periodic
// inflow boundary (PCPhysBCType encoding: inflow == 1); dst is the would-be
// destination (outside the domain for a boundary face)
AMREX_FORCE_INLINE
bool
is_inflow_boundary_face(
  const amrex::IntVect& dst,
  const amrex::Box& domain,
  const amrex::IntVect& is_periodic,
  const amrex::IntVect& lo_bc,
  const amrex::IntVect& hi_bc)
{
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    if (is_periodic[dir] != 0) {
      continue;
    }
    if (dst[dir] < domain.smallEnd(dir) && lo_bc[dir] == 1) {
      return true;
    }
    if (dst[dir] > domain.bigEnd(dir) && hi_bc[dir] == 1) {
      return true;
    }
  }
  return false;
}

} // namespace boundary
} // namespace clem
#endif
