#ifndef CLEM_ADVECTION_H
#define CLEM_ADVECTION_H

#include <AMReX.H>


namespace clem{
    namespace advection {
        enum class Face {Left=0, Right=1, Bottom=2, Top=3, South=4, North=5};

        struct FluxDir {
            amrex::Real flux;
            Face face;
        };

        AMREX_GPU_DEVICE
        constexpr amrex::IntVect face_offsets[] = {
            amrex::IntVect{AMREX_D_DECL(-1, 0, 0)}, // Left
            amrex::IntVect{AMREX_D_DECL(1, 0, 0)}, // Right
            amrex::IntVect{AMREX_D_DECL(0,-1, 0)}, // Bottom
            amrex::IntVect{AMREX_D_DECL(0, 1, 0)}, // Top
        #if AMREX_SPACEDIM == 3
            amrex::IntVect{AMREX_D_DECL(0, 0,-1)}, // South
            amrex::IntVect{AMREX_D_DECL(0, 0, 1)}, // North
        #endif
        };

        AMREX_GPU_DEVICE
        inline
        amrex::IntVect neighbor_iv(const amrex::IntVect& iv, Face face)
        {
            return iv + face_offsets[static_cast<int>(face)];
        }

    }

    namespace boundaryCondition{
        //To apply periodic BC - wraps cell indices based on periodic geometry
        AMREX_FORCE_INLINE
        std::pair<amrex::IntVect, bool> wrap_iv(const amrex::IntVect& iv, 
                                                const amrex::Box& domain,
                                                const amrex::IntVect& is_periodic)
        {
            amrex::IntVect r = iv;
            bool was_wrapped = false;
            
            bool is_inside_the_box = domain.contains(r);
            if (is_inside_the_box) {
                return {r, false};
            } else {
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    if (is_periodic[dir]) {
                        int length = domain.length(dir);
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
        }
        
        //outflow implementation - check if cell is outside domain on a non-periodic boundary
        AMREX_FORCE_INLINE
        bool is_outside_domain_on_outflow(const amrex::IntVect& iv,
                                          const amrex::Box& domain,
                                          const amrex::IntVect& is_periodic,
                                          const amrex::IntVect& lo_bc,
                                          const amrex::IntVect& hi_bc)
        {
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                if (is_periodic[dir]) continue;  // Skip periodic directions
                
                // outflow implementation - check lower boundary
                if (iv[dir] < domain.smallEnd(dir)) {
                    // Check if this direction has outflow BC at low boundary
                    // Note: BC type 2 is outflow in PeleC (PCPhysBCType::outflow)
                    if (lo_bc[dir] == 2) return true;  // Particle escaped on outflow
                }
                
                // outflow implementation - check upper boundary
                if (iv[dir] > domain.bigEnd(dir)) {
                    // Check if this direction has outflow BC at high boundary
                    if (hi_bc[dir] == 2) return true;  // Particle escaped on outflow
                }
            }
            return false;
        }

        //inflow implementation - check if crossing a face enters the domain through a non-periodic inflow boundary.
        //dst is the would-be destination cell (outside the domain for a boundary face). BC type 1 is inflow in PeleC (PCPhysBCType::inflow).
        AMREX_FORCE_INLINE
        bool is_inflow_boundary_face(const amrex::IntVect& dst,
                                     const amrex::Box& domain,
                                     const amrex::IntVect& is_periodic,
                                     const amrex::IntVect& lo_bc,
                                     const amrex::IntVect& hi_bc)
        {
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                if (is_periodic[dir]) continue;  // periodic directions wrap, they do not inject mass

                if (dst[dir] < domain.smallEnd(dir) && lo_bc[dir] == 1) return true;  // inflow at low boundary
                if (dst[dir] > domain.bigEnd(dir)   && hi_bc[dir] == 1) return true;  // inflow at high boundary
            }
            return false;
        }
    }
}

#endif  /*CLEM_ADVECTION_H*/