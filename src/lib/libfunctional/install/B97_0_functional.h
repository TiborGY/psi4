#ifndef B97_0_functional_h
#define B97_0_functional_h
/**********************************************************
* B97_0_functional.h: declarations for B97_0_functional for KS-DFT
* Robert Parrish, robparrish@gmail.com
* Autogenerated by MATLAB Script on 06-Mar-2012
*
***********************************************************/
#include "functional.h"

namespace psi { namespace functional {

class B97_0_Functional : public Functional {
public:
    B97_0_Functional(int npoints, int deriv);
    virtual ~B97_0_Functional();
    virtual void computeRKSFunctional(boost::shared_ptr<RKSFunctions> prop);
    virtual void computeUKSFunctional(boost::shared_ptr<RKSFunctions> prop);
};
}}
#endif

