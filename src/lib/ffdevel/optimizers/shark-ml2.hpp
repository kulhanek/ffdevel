#ifndef SHARK_ML2_HEADER
#define SHARK_ML2_HEADER

#include <shark/ObjectiveFunctions/AbstractObjectiveFunction.h>
#include <shark/Core/Random.h>

// this is for NB2LJ conversion

// optimized function ----------------------------

struct COptSharkFce2 : public shark::SingleObjectiveFunction {
public:
    COptSharkFce2(void);
    std::string name(void) const;
    std::size_t numberOfVariables(void) const;
    double eval(SearchPointType const& x) const;
};

#endif
