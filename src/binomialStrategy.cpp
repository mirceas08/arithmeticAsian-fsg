#ifndef BINOMIALSTRATEGY_CPP
#define BINOMIALSTRATEGY_CPP

#include "binomialStrategy.h"
#include <cmath>


/* ------------------------ Base class ------------------------ */

BinomialStrategy::BinomialStrategy(double _sigma, double _r, double _dt)
{
    sigma = _sigma;
    r = _r;
    dt = _dt;
}

BinomialStrategy::~BinomialStrategy() {}


/* ------------------------ Cox-Ross-Rubinstein  ------------------------ */

CRR::CRR(double sigma, double r, double dt):
    BinomialStrategy(sigma, r, dt)
{
    u = std::exp(sigma * std::sqrt(dt));
    d = std::exp(-sigma * std::sqrt(dt));
    p = (std::exp(r*dt) - d) / (u-d);
}

CRR::~CRR() {}


/* ------------------------ Jarrow-Rudd  ------------------------ */

JarrowRudd::JarrowRudd(double sigma, double r, double dt):
    BinomialStrategy(sigma, r, dt)
{
    double drift = r - 0.5 * sigma * sigma;
    u = std::exp(drift*dt + sigma * std::sqrt(dt));
    d = std::exp(drift*dt - sigma * std::sqrt(dt));
    p = 0.5;
}

JarrowRudd::~JarrowRudd() {}

#endif // BINOMIALSTRATEGY_CPP
