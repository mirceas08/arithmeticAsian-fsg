#ifndef BINOMIALSTRATEGY_CPP
#define BINOMIALSTRATEGY_CPP

#include "binomialStrategy.h"
#include <cmath>

/* ------------------------------------------------------------ */
/* ------------------------ Base class ------------------------ */
/* ------------------------------------------------------------ */

BinomialStrategy::BinomialStrategy(double _sigma, double _r, double _dt)
{
    sigma = _sigma;
    r = _r;
    dt = _dt;
}

BinomialStrategy::~BinomialStrategy() {}

/* ----------------------------------------------------------------- */
/* ------------------------ Derived classes ------------------------ */
/* ----------------------------------------------------------------- */


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


/* ------------------------ Tian  ------------------------ */

Tian::Tian(double sigma, double r, double dt):
    BinomialStrategy(sigma, r, dt)
{
    double v = std::exp(sigma*sigma*dt);
    double a = std::exp(r*dt);
    double sqrtV = std::sqrt(v*v + 2*v - 3);

    u = 0.5 * a * v * (v + 1 + sqrtV);
    d = 0.5 * a * v * (v + 1 - sqrtV);
    p = (a-d) / (u-d);
}

Tian::~Tian() {}

/* ------------------------ Jarrow-Rudd Risk Neutral ------------------------ */

JarrowRudd_RN::JarrowRudd_RN(double sigma, double r, double dt):
    BinomialStrategy(sigma, r, dt)
{
    double drift = r - 0.5 * sigma * sigma;
    u = std::exp(drift*dt + sigma * std::sqrt(dt));
    d = std::exp(drift*dt - sigma * std::sqrt(dt));
    p = (std::exp(r*dt) - d) / (u-d);
}

JarrowRudd_RN::~JarrowRudd_RN() {}


/* ------------------------ Cox-Ross-Rubinstein with drift ------------------------ */

CRR_drift::CRR_drift(double sigma, double r, double dt, double _eta):
    BinomialStrategy(sigma, r, dt), eta(_eta)
{
    u = std::exp(eta*dt + sigma * std::sqrt(dt));
    d = std::exp(eta*dt -sigma * std::sqrt(dt));
    p = (std::exp(r*dt) - d) / (u-d);
}

CRR_drift::~CRR_drift() {}

/* ------------------------ Leisen-Reimer ------------------------ */

double h(double x, int n);

LeisenReimer::LeisenReimer(double sigma, double r, double dt, double _S0, double _K, double _T, int _n):
    BinomialStrategy(sigma, r, dt), S0(_S0), K(_K), T(_T), n(_n)
{
    double sqrtT = sigma * std::sqrt(T);
    double d1 = (std::log(S0/K) + (r + 0.5 * sigma * sigma) * T) / sqrtT;
    double d2 = d1 - sqrtT;

    p_bar = h(d1, n);
    p = h(d2, n);

    u = std::exp(r*dt) * (p_bar / p);
    d = (std::exp(r*dt) - p * u) / (1.0 - p);
}

LeisenReimer::~LeisenReimer() {}


#endif // BINOMIALSTRATEGY_CPP
