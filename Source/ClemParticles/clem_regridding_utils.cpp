#include "ClemParticles/clem_regridding_utils.h"


namespace clem{
    
    /**
     * @brief Reset the accumulator data structure to zero. 
     * This is useful to initialize the structure before accumulating new data.
     */
    void regridding::LemAccumulator::reset()
    {
        mass     = 0.0;
        pressure = 0.0;
        volume   = 0.0;
        density  = 0.0;
        energy   = 0.0;
        temp     = 0.0;

        for (auto& v : sp) {
            v = 0.0;
        }
    }

    void  regridding::AccumulateDataIntoStructureGivenMassPortion(LemAccumulator& acc, const ClemParticles::ParticleType& p, amrex::Real mass_need)
        {
            acc.mass      += mass_need;
            acc.pressure  += p.rdata(RealData::pressure) * mass_need;
            acc.volume    += mass_need / p.rdata(RealData::density);
            //internal energy updated based on temperrature
            //amrex::Real energy    = 0.0;
            acc.temp      += p.rdata(RealData::temp) * mass_need;

            amrex::Real one_rho = 1.0 / p.rdata(RealData::density);

            //m*Y .. Here is transformed into primitves
            for(int sp = 0; sp < NUM_SPECIES ; sp++)
                acc.sp[sp] += mass_need * p.rdata(RealData::SP(sp)) * one_rho;
        }

    void  regridding::AccumulateDataIntoStructureGivenMassPortion(LemAccumulator& acc, const LemAccumulator& aux, amrex::Real mass_need)
        {
            acc.mass      += mass_need;
            acc.pressure  += aux.pressure * mass_need;
            acc.volume    += mass_need / aux.density;
            //internal energy updated based on temperrature
            //amrex::Real energy    = 0.0;
            acc.temp      += aux.temp * mass_need;

            amrex::Real one_rho = 1.0 / aux.density;

            //m*Y .. Here is transformed into primitves
            for(int sp = 0; sp < NUM_SPECIES ; sp++)
                acc.sp[sp] += mass_need * aux.sp[sp] * one_rho;
        }

    void regridding::StoreRemainningDataIntoAux(LemAccumulator& data, const ClemParticles::ParticleType& p )
        {   
            //Store only the relevant variables
            data.mass       = p.rdata(RealData::mass);
            data.pressure   = p.rdata(RealData::pressure);
            data.density    = p.rdata(RealData::density);
            data.volume     = p.rdata(RealData::volume);
            data.temp       = p.rdata(RealData::temp);

            //Keep conservatives
            for(int sp = 0; sp < NUM_SPECIES; sp++)
                data.sp[sp] = p.rdata(RealData::SP(sp));
        }

    void regridding::UpdateLemFromAccumulator(const LemAccumulator& data, ClemParticles::ParticleType& dst, amrex::Real new_dm){

            dst.rdata(RealData::mass)           = data.mass ;
            dst.rdata(RealData::pressure)       = data.pressure / new_dm    ;
            dst.rdata(RealData::volume)         = data.volume   ;
            dst.rdata(RealData::density)        = data.mass / data.volume   ;
            dst.rdata(RealData::temp)           = data.temp /new_dm ;

            amrex::Real one_vol                 = 1.0 / data.volume;

            //data contains m*Y then to tranlsate to conservative m*y/vol
            for(int sp = 0; sp < NUM_SPECIES; sp++)
                dst.rdata(RealData::SP(sp))     = data.sp[sp] * one_vol;
        }

    void regridding::AccumulateDataIntoStructure(LemAccumulator& acc, LemAccumulator& aux){
            AccumulateDataIntoStructureGivenMassPortion(acc, aux, aux.mass);
        }

    void regridding::AccumulateDataIntoStructure(LemAccumulator& acc, ClemParticles::ParticleType& p){
            AccumulateDataIntoStructureGivenMassPortion(acc, p, p.rdata(RealData::mass));
        }



}