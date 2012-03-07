#ifndef EDF1_C_functional_h
#define EDF1_C_functional_h
/**********************************************************
* EDF1_C_functional.h: declarations for EDF1_C_functional for KS-DFT
* Robert Parrish, robparrish@gmail.com
* Autogenerated by MATLAB Script on 06-Mar-2012
*
***********************************************************/
#include "functional.h"

namespace psi { namespace functional {

class EDF1_C_Functional : public Functional {
public:
    EDF1_C_Functional(int npoints, int deriv);
    virtual ~EDF1_C_Functional();
    virtual void computeRKSFunctional(boost::shared_ptr<RKSFunctions> prop);
    virtual void computeUKSFunctional(boost::shared_ptr<RKSFunctions> prop);
};
}}
#endif

