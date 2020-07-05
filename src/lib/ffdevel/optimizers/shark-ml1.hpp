#ifndef SHARK_ML1_HEADER
#define SHARK_ML1_HEADER

#include <shark/ObjectiveFunctions/AbstractObjectiveFunction.h>
#include <shark/Core/Random.h>

// this is for FF parameter optimization

// methods ---------------------------------------
const int           SHARK_CMA_ES = 1;

// optimized function ----------------------------

struct COptSharkFce1 : public shark::SingleObjectiveFunction {
public:
    COptSharkFce1(void);
    std::string name(void) const;
    std::size_t numberOfVariables(void) const;
    double eval(SearchPointType const& x) const;
};

#endif
