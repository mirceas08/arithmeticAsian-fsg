#ifndef OPTION_CPP
#define OPTION_CPP

#include "option.h"

/* ------------------------ Option base class ------------------------ */

Option::Option() {}

Option::~Option() {}

/* ------------------------ Average call option ------------------------ */

AverageCallOption::AverageCallOption(vec _K):
    K(_K) {}

AverageCallOption::~AverageCallOption() {}

vec AverageCallOption::payoff(const vec &S_average)
{
    int vec_size = S_average.size();
    return max(S_average - K, zeros<vec>(vec_size));
}

/* ------------------------ Average put option ------------------------ */

AveragePutOption::AveragePutOption(vec _K):
    K(_K) {}

AveragePutOption::~AveragePutOption() {}

vec AveragePutOption::payoff(const vec &S_average)
{
    int vec_size = S_average.size();
    return max(K - S_average, zeros<vec>(vec_size));
}

#endif // OPTION_CPP
