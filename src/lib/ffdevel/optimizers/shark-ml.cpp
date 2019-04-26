#include <shark-ml.hpp>
#include <memory>
#include <chrono>
#include <shark/Algorithms/DirectSearch/CMA.h>

using namespace shark;
using namespace std;

// -----------------------------------------------------------------------------

int                 OptMethod;
int                 NumOfParams;
shared_ptr<CMA>     CMAOpt;
COptSharkFce        OptSharkFce;
shared_ptr<double>  tmp_x;
std::mt19937        RnG;

// -----------------------------------------------------------------------------

extern "C" void shark_create_(int* nactparms,int* method,double* initial_step,int* rngseed,double* initial_params)
{
    NumOfParams = *nactparms;
    RealVector params(NumOfParams);
    for(size_t i=0; i < NumOfParams; i++){
            params(i) = initial_params[i];
    }
    tmp_x = shared_ptr<double>(new double[NumOfParams]);
    if( *rngseed == 0 ) {
        unsigned seed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch()).count();
        RnG.seed(seed);
    } else {
        // setup RNG
        RnG.seed(*rngseed);
    }

    // currently only CMA is supported
    OptMethod = *method;
    switch(OptMethod){
        case SHARK_CMA_ES:
            CMAOpt = shared_ptr<CMA>(new CMA(RnG));
            CMAOpt->setInitialSigma(*initial_step);
            CMAOpt->init(OptSharkFce,params);
            break;
    }
}

// -----------------------------------------------------------------------------

extern "C" void shark_dostep_(double* error)
{
    switch(OptMethod){
        case SHARK_CMA_ES:
            CMAOpt->step(OptSharkFce);
            *error = CMAOpt->solution().value;
            break;
    }
}

// -----------------------------------------------------------------------------

extern "C" void shark_getsol_(double* params)
{
    switch(OptMethod){
        case SHARK_CMA_ES:
            for(size_t i=0; i < NumOfParams; i++ ){
                params[i] = CMAOpt->solution().point(i);
            }
            break;
    }
}

// -----------------------------------------------------------------------------

extern "C" void shark_destroy_(void)
{
    // destroy data
    switch(OptMethod){
        case SHARK_CMA_ES:
            CMAOpt = shared_ptr<CMA>();
            break;
    }
    tmp_x = shared_ptr<double>();
}

// -----------------------------------------------------------------------------

COptSharkFce::COptSharkFce(void)
{
    m_features |= HAS_VALUE;
}

// -----------------------------------------------------------------------------

std::string COptSharkFce::name(void) const
{
    return("FFDevel");
}

// -----------------------------------------------------------------------------

std::size_t COptSharkFce::numberOfVariables(void) const
{
    return(NumOfParams);
}

// -----------------------------------------------------------------------------

extern "C" void opt_shark_fce_(int* nprms,double* prms,double* value);

double COptSharkFce::eval(SearchPointType const& x) const
{
    double* prms = tmp_x.get();
    for(size_t i=0; i < NumOfParams; i++) prms[i] = x(i);
    double  value = 0.0;

    // call fortran method
    opt_shark_fce_(&NumOfParams,prms,&value);

    return(value);
}

// -----------------------------------------------------------------------------
