#ifndef OPTION_H
#define OPTION_H

#include <armadillo>
using namespace arma;

// base class
class Option
{
public:
    Option();
    virtual ~Option();

    virtual vec payoff(const vec &S_average) = 0;
};

// derived classes

/* ------------------------ Average options ------------------------ */

class AverageCallOption: public Option
{
public:
    vec K;
public:
    AverageCallOption(vec _K);
    virtual ~AverageCallOption();

    virtual vec payoff(const vec &S_average);
};

class AveragePutOption: public Option
{
public:
    vec K;
public:
    AveragePutOption(vec _K);
    virtual ~AveragePutOption();

    virtual vec payoff(const vec &S_average);
};

#endif // OPTION_H
