#include <shark-ml2.hpp>
#include <memory>
#include <chrono>
#include <shark/Algorithms/DirectSearch/CMA.h>

using namespace shark;
using namespace std;

// -----------------------------------------------------------------------------

int                 NumOfParams2;
shared_ptr<CMA>     CMAOpt2;
COptSharkFce2       OptSharkFce2;
shared_ptr<double>  tmp_x2;
std::mt19937        RnG2;

// -----------------------------------------------------------------------------

extern "C" void shark_create2_(int* nactparms,double* initial_step,double* initial_params)
{
    NumOfParams2 = *nactparms;
    if( NumOfParams2 <= 0 ) return;

    RealVector params(NumOfParams2);
    for(size_t i=0; i < NumOfParams2; i++){
            params(i) = initial_params[i];
    }
    tmp_x2 = shared_ptr<double>(new double[NumOfParams2],std::default_delete<double[]>());

    // currently only CMA is supported
    CMAOpt2 = shared_ptr<CMA>(new CMA(RnG2));
    CMAOpt2->setInitialSigma(*initial_step);
    CMAOpt2->init(OptSharkFce2,params);
}

// -----------------------------------------------------------------------------

extern "C" void shark_set_rngseed2_(int* seed)
{
    if( *seed == 0 ) {
        unsigned seed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch()).count();
        RnG2.seed(seed);
    } else {
        // setup RNG
        RnG2.seed(*seed);
    }
}

// -----------------------------------------------------------------------------

extern "C" void shark_dostep2_(double* error)
{
    *error = 0.0;
    if( NumOfParams2 <= 0 ) return;
    CMAOpt2->step(OptSharkFce2);
    *error = CMAOpt2->solution().value;
}

// -----------------------------------------------------------------------------

extern "C" void shark_getsol2_(double* params)
{
    if( NumOfParams2 <= 0 ) return;
    for(size_t i=0; i < NumOfParams2; i++ ){
        params[i] = CMAOpt2->solution().point(i);
    }
}

// -----------------------------------------------------------------------------

extern "C" void shark_destroy2_(void)
{
    // destroy data
    CMAOpt2 = shared_ptr<CMA>();
    tmp_x2 = shared_ptr<double>();
}

// -----------------------------------------------------------------------------

COptSharkFce2::COptSharkFce2(void)
{
    m_features |= HAS_VALUE;
}

// -----------------------------------------------------------------------------

std::string COptSharkFce2::name(void) const
{
    return("FFDevel");
}

// -----------------------------------------------------------------------------

std::size_t COptSharkFce2::numberOfVariables(void) const
{
    return(NumOfParams2);
}

// -----------------------------------------------------------------------------

extern "C" void opt_shark_fce2_(int* nprms,double* prms,double* value);

double COptSharkFce2::eval(SearchPointType const& x) const
{
    double* prms = tmp_x2.get();
    for(size_t i=0; i < NumOfParams2; i++) prms[i] = x(i);
    double  value = 0.0;

    // call fortran method
    opt_shark_fce2_(&NumOfParams2,prms,&value);

    return(value);
}

// -----------------------------------------------------------------------------
