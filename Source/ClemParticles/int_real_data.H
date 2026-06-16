#ifndef INT_REAL_DATA_H
#define INT_REAL_DATA_H

#include "ClemParticles/clem_index_definition.H"

namespace clem{
//==================================================
// Particle real data layout
//==================================================
    struct RealData
    {
        enum{
            mass = 0, //mass
            pressure, //1 pressure
            volume, //2 volume
            flux, //3// Cosnervatives variables .. Flux
            density, //4 density
            energy, //rho_eint 5 ueint
            temp, //6  //Species  temperature
            sp_base, 
            ncomps = sp_base + NUM_SPECIES
            };

            static constexpr int SP(int n) {
            return sp_base + n;
            }
    };
//==================================================
// Particle integer data layout
//==================================================
    struct IntData
    {       
        enum{ idx = 0, //Relative position within the LemArray 
                i,     // cell idx (i,j ,k)     i
                j,      //j
                k,      //k
                advected, // is advected
                ncomps};
    };
}
#endif /* INT_REAL_DATA_H */