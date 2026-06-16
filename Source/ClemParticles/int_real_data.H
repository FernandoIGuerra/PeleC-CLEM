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
            mass = 0,
            pressure, //1
            volume, //2
            flux, //3// Cosnervatives variables
            density, //4
            energy, //rho_eint 5
            temp, //6  //Species 
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
                i,     // cell idx (i,j ,k)     
                j, 
                k,
                advected, 
                ncomps};
    };
}
#endif /* INT_REAL_DATA_H */