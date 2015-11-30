#ifndef OPTION_CPP
#define OPTION_CPP

#include "option.h"

/* ------------------------ Option base class ------------------------ */

Option::Option() {}

Option::~Option() {}

/* ------------------------ Average call option ------------------------ */

AverageCallOption::AverageCallOption() {}

AverageCallOption::~AverageCallOption() {}

vec AverageCallOption::payoff(const vec &S, const vec &strike)
{
    int vec_size = S.size();
    return max(S - strike, zeros<vec>(vec_size));
}

/* ------------------------ Average put option ------------------------ */

AveragePutOption::AveragePutOption() {}

AveragePutOption::~AveragePutOption() {}

vec AveragePutOption::payoff(const vec &S, const vec &strike)
{
    int vec_size = S.size();
    return max(strike - S, zeros<vec>(vec_size));
}

/* ####################################################################### */

/* ------------------------ Asian call option ------------------------ */

AsianCallOption::AsianCallOption() {}

AsianCallOption::~AsianCallOption() {}

vec AsianCallOption::payoff(const vec &S, const vec &strike)
{
    int vec_size = S.size();
    return max(S - strike, zeros<vec>(vec_size));
}

/* ------------------------ Average put option ------------------------ */

AsianPutOption::AsianPutOption() {}

AsianPutOption::~AsianPutOption() {}

vec AsianPutOption::payoff(const vec &S, const vec &strike)
{
    int vec_size = S.size();
    return max(strike - S, zeros<vec>(vec_size));
}

#endif // OPTION_CPP
