#ifndef BINOMIALSTRATEGY_H
#define BINOMIALSTRATEGY_H

// base class
class BinomialStrategy
{
public:
    double u;
    double d;
    double p;

    double r;
    double sigma;
    double dt;
public:
    BinomialStrategy(double sigma, double r, double dt);
    virtual ~BinomialStrategy();
};

// Cox-Ross-Rubinstein
class CRR: public BinomialStrategy
{
public:
    CRR(double sigma, double r, double dt);
    virtual ~CRR();
};

// Jarrow-Rudd
class JarrowRudd: public BinomialStrategy
{
public:
    JarrowRudd(double sigma, double r, double dt);
    virtual ~JarrowRudd();
};

#endif // BINOMIALSTRATEGY_H
