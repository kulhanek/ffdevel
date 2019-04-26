#ifndef SHARK_ML_HEADER
#define SHARK_ML_HEADER

#include <shark/ObjectiveFunctions/AbstractObjectiveFunction.h>
#include <shark/Core/Random.h>

// methods ---------------------------------------
const int           SHARK_CMA_ES = 1;

// optimized function ----------------------------

struct COptSharkFce : public shark::SingleObjectiveFunction {
public:
    COptSharkFce(void);
    std::string name(void) const;
    std::size_t numberOfVariables(void) const;
    double eval(SearchPointType const& x) const;
};

#endif
