#include <shark-ml1.hpp>
#include <memory>
#include <chrono>
#include <shark/Algorithms/DirectSearch/CMA.h>

using namespace shark;
using namespace std;

// -----------------------------------------------------------------------------

int                 OptMethod1;
int                 NumOfParams1;
shared_ptr<CMA>     CMAOpt1;
COptSharkFce1       OptSharkFce1;
shared_ptr<double>  tmp_x1;
std::mt19937        RnG1;

// -----------------------------------------------------------------------------

extern "C" void shark_create1_(int* nactparms,int* method,double* initial_step,double* initial_params)
{
    NumOfParams1 = *nactparms;
    if( NumOfParams1 <= 0 ) return;

    RealVector params(NumOfParams1);
    for(size_t i=0; i < NumOfParams1; i++){
            params(i) = initial_params[i];
    }
    tmp_x1 = shared_ptr<double>(new double[NumOfParams1]);

    // currently only CMA is supported
    OptMethod1 = *method;
    switch(OptMethod1){
        case SHARK_CMA_ES:
            CMAOpt1 = shared_ptr<CMA>(new CMA(RnG1));
            CMAOpt1->setInitialSigma(*initial_step);
            CMAOpt1->init(OptSharkFce1,params);
            break;
    }
}

// -----------------------------------------------------------------------------

extern "C" void shark_set_rngseed1_(int* seed)
{
    if( *seed == 0 ) {
        unsigned seed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch()).count();
        RnG1.seed(seed);
    } else {
        // setup RNG
        RnG1.seed(*seed);
    }
}

// -----------------------------------------------------------------------------

extern "C" void shark_dostep1_(double* error)
{
    *error = 0.0;
    if( NumOfParams1 <= 0 ) return;
    switch(OptMethod1){
        case SHARK_CMA_ES:
            CMAOpt1->step(OptSharkFce1);
            *error = CMAOpt1->solution().value;
            break;
    }
}

// -----------------------------------------------------------------------------

extern "C" void shark_getsol1_(double* params)
{
    if( NumOfParams1 <= 0 ) return;
    switch(OptMethod1){
        case SHARK_CMA_ES:
            for(size_t i=0; i < NumOfParams1; i++ ){
                params[i] = CMAOpt1->solution().point(i);
            }
            break;
    }
}

// -----------------------------------------------------------------------------

extern "C" void shark_destroy1_(void)
{
    // destroy data
    switch(OptMethod1){
        case SHARK_CMA_ES:
            CMAOpt1 = shared_ptr<CMA>();
            break;
    }
    tmp_x1 = shared_ptr<double>();
}

// -----------------------------------------------------------------------------

COptSharkFce1::COptSharkFce1(void)
{
    m_features |= HAS_VALUE;
}

// -----------------------------------------------------------------------------

std::string COptSharkFce1::name(void) const
{
    return("FFDevel");
}

// -----------------------------------------------------------------------------

std::size_t COptSharkFce1::numberOfVariables(void) const
{
    return(NumOfParams1);
}

// -----------------------------------------------------------------------------

extern "C" void opt_shark_fce1_(int* nprms,double* prms,double* value);

double COptSharkFce1::eval(SearchPointType const& x) const
{
    double* prms = tmp_x1.get();
    for(size_t i=0; i < NumOfParams1; i++) prms[i] = x(i);
    double  value = 0.0;

    // call fortran method
    opt_shark_fce1_(&NumOfParams1,prms,&value);

    return(value);
}

// -----------------------------------------------------------------------------
