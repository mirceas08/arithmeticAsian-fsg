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

    virtual vec payoff(const vec &myVector, const double &myDouble) = 0;
};

// derived classes

/* ------------------------ Average options ------------------------ */

class AverageCallOption: public Option
{
public:
    AverageCallOption();
    virtual ~AverageCallOption();

    virtual vec payoff(const vec &myVector, const double &myDouble);
};

class AveragePutOption: public Option
{
public:
    AveragePutOption();
    virtual ~AveragePutOption();

    virtual vec payoff(const vec &myVector, const double &myDouble);
};

/* ------------------------ Asian options ------------------------ */

class AsianCallOption: public Option
{
public:
    AsianCallOption();
    virtual ~AsianCallOption();

    virtual vec payoff(const vec &myVector, const double &myDouble);
};

class AsianPutOption: public Option
{
public:
    AsianPutOption();
    virtual ~AsianPutOption();

    virtual vec payoff(const vec &myVector, const double &myDouble);
};

#endif // OPTION_H
