#ifndef OPTION_CPP
#define OPTION_CPP

#include "option.h"

/* ------------------------ Option base class ------------------------ */

Option::Option() {}

Option::~Option() {}

/* ------------------------ Average call option ------------------------ */

AverageCallOption::AverageCallOption() {}

AverageCallOption::~AverageCallOption() {}

vec AverageCallOption::payoff(const vec &myVector, const double &myDouble)
{
    int vecSize = myVector.size();
    return max(myVector - myDouble, zeros<vec>(vecSize));
}

/* ------------------------ Average put option ------------------------ */

AveragePutOption::AveragePutOption() {}

AveragePutOption::~AveragePutOption() {}

vec AveragePutOption::payoff(const vec &myVector, const double &myDouble)
{
    int vecSize = myVector.size();
    return max(myDouble - myVector, zeros<vec>(vecSize));
}

/* ####################################################################### */

/* ------------------------ Asian call option ------------------------ */

AsianCallOption::AsianCallOption() {}

AsianCallOption::~AsianCallOption() {}

vec AsianCallOption::payoff(const vec &myVector, const double &myDouble)
{
    int vecSize = myVector.size();
    return max(myDouble - myVector, zeros<vec>(vecSize));
}

/* ------------------------ Average put option ------------------------ */

AsianPutOption::AsianPutOption() {}

AsianPutOption::~AsianPutOption() {}

vec AsianPutOption::payoff(const vec &myVector, const double &myDouble)
{
    int vecSize = myVector.size();
    return max(myVector - myDouble, zeros<vec>(vecSize));
}

#endif // OPTION_CPP
