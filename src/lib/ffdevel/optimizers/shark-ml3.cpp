#include <shark-ml3.hpp>
#include <memory>
#include <chrono>
#include <shark/Algorithms/DirectSearch/CMA.h>

using namespace shark;
using namespace std;

// -----------------------------------------------------------------------------

int                 NumOfParams3;
shared_ptr<CMA>     CMAOpt3;
COptSharkFce3       OptSharkFce3;
shared_ptr<double>  tmp_x3;
std::mt19937        RnG3;

// -----------------------------------------------------------------------------

extern "C" void shark_create3_(int* nactparms,double* initial_step,double* initial_params)
{
    NumOfParams3 = *nactparms;
    if( NumOfParams3 <= 0 ) return;

    RealVector params(NumOfParams3);
    for(size_t i=0; i < NumOfParams3; i++){
            params(i) = initial_params[i];
    }
    tmp_x3 = shared_ptr<double>(new double[NumOfParams3]);

    // currently only CMA is supported
    CMAOpt3 = shared_ptr<CMA>(new CMA(RnG3));
    CMAOpt3->setInitialSigma(*initial_step);
    CMAOpt3->init(OptSharkFce3,params);
}

// -----------------------------------------------------------------------------

extern "C" void shark_set_rngseed3_(int* seed)
{
    if( *seed == 0 ) {
        unsigned seed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch()).count();
        RnG3.seed(seed);
    } else {
        // setup RNG
        RnG3.seed(*seed);
    }
}

// -----------------------------------------------------------------------------

extern "C" void shark_dostep3_(double* error)
{
    *error = 0.0;
    if( NumOfParams3 <= 0 ) return;
    CMAOpt3->step(OptSharkFce3);
    *error = CMAOpt3->solution().value;
}

// -----------------------------------------------------------------------------

extern "C" void shark_getsol3_(double* params)
{
    if( NumOfParams3 <= 0 ) return;
    for(size_t i=0; i < NumOfParams3; i++ ){
        params[i] = CMAOpt3->solution().point(i);
    }
}

// -----------------------------------------------------------------------------

extern "C" void shark_destroy3_(void)
{
    // destroy data
    CMAOpt3 = shared_ptr<CMA>();
    tmp_x3 = shared_ptr<double>();
}

// -----------------------------------------------------------------------------

COptSharkFce3::COptSharkFce3(void)
{
    m_features |= HAS_VALUE;
}

// -----------------------------------------------------------------------------

std::string COptSharkFce3::name(void) const
{
    return("FFDevel");
}

// -----------------------------------------------------------------------------

std::size_t COptSharkFce3::numberOfVariables(void) const
{
    return(NumOfParams3);
}

// -----------------------------------------------------------------------------

extern "C" void opt_shark_fce3_(int* nprms,double* prms,double* value);

double COptSharkFce3::eval(SearchPointType const& x) const
{
    double* prms = tmp_x3.get();
    for(size_t i=0; i < NumOfParams3; i++) prms[i] = x(i);
    double  value = 0.0;

    // call fortran method
    opt_shark_fce3_(&NumOfParams3,prms,&value);

    return(value);
}

// -----------------------------------------------------------------------------
