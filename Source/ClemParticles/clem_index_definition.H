#ifndef CLEM_INDEX
#define CLEM_INDEX

#include "PelePhysics.H"
#include "mechanism.H"

#define NUM_LEM 9
#define CLEM_URHO 0
#define CLEM_UTEMP 1
#define CLEM_UEINT 2
#define CLEM_UFS 3
#define DO_CLEM_PARTICLE 1

#define CLEM_UPRES 1
#define CLEM_NVAR (CLEM_UFS + NUM_SPECIES)
//for trasnport properties

#define CLEM_dComp_lambda (NUM_SPECIES)


namespace clem{
    namespace constants{
        constexpr amrex::Real isentropicGamma               = 1.4;
        constexpr amrex::Real isentropicRatio               = 0.4/1.4; 
    }
}
#endif