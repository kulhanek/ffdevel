#ifndef SHARK_ML3_HEADER
#define SHARK_ML3_HEADER

#include <shark/ObjectiveFunctions/AbstractObjectiveFunction.h>
#include <shark/Core/Random.h>

// this is for QNB isovalue calculations

// optimized function ----------------------------

struct COptSharkFce3 : public shark::SingleObjectiveFunction {
public:
    COptSharkFce3(void);
    std::string name(void) const;
    std::size_t numberOfVariables(void) const;
    double eval(SearchPointType const& x) const;
};

#endif
