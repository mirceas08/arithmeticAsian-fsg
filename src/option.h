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

    virtual vec payoff(const vec &S, const vec &strike) = 0;
};

// derived classes

/* ------------------------ Average options ------------------------ */

class AverageCallOption: public Option
{
public:
    AverageCallOption();
    virtual ~AverageCallOption();

    virtual vec payoff(const vec &S, const vec &strike);
};

class AveragePutOption: public Option
{
public:
    AveragePutOption();
    virtual ~AveragePutOption();

    virtual vec payoff(const vec &S, const vec &strike);
};

/* ------------------------ Asian options ------------------------ */

class AsianCallOption: public Option
{
public:
    AsianCallOption();
    virtual ~AsianCallOption();

    virtual vec payoff(const vec &S, const vec &strike);
};

class AsianPutOption: public Option
{
public:
    AsianPutOption();
    virtual ~AsianPutOption();

    virtual vec payoff(const vec &S, const vec &strike);
};

#endif // OPTION_H
