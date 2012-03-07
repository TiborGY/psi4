#ifndef wB_X_functional_h
#define wB_X_functional_h
/**********************************************************
* wB_X_functional.h: declarations for wB_X_functional for KS-DFT
* Robert Parrish, robparrish@gmail.com
* Autogenerated by MATLAB Script on 07-Mar-2012
*
***********************************************************/
#include "functional.h"

namespace psi { namespace functional {

class wB_X_Functional : public Functional {
public:
    wB_X_Functional(int npoints, int deriv);
    virtual ~wB_X_Functional();
    virtual void computeRKSFunctional(boost::shared_ptr<RKSFunctions> prop);
    virtual void computeUKSFunctional(boost::shared_ptr<RKSFunctions> prop);
};
}}
#endif

