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
    BinomialStrategy(double _sigma, double _r, double _dt);
    virtual ~BinomialStrategy();
};

// Cox-Ross-Rubinstein (CRR)
class CRR: public BinomialStrategy
{
public:
    CRR(double sigma, double r, double dt);
    virtual ~CRR();
};

// Jarrow-Rudd (JR)
class JarrowRudd: public BinomialStrategy
{
public:
    JarrowRudd(double sigma, double r, double dt);
    virtual ~JarrowRudd();
};

// Tian (TIAN)
class Tian: public BinomialStrategy
{
public:
    Tian(double sigma, double r, double dt);
    virtual ~Tian();
};

// Jarrow-Rudd Risk Neutral (JR-RN)
class JarrowRudd_RN: public BinomialStrategy
{
public:
    JarrowRudd_RN(double sigma, double r, double dt);
    virtual ~JarrowRudd_RN();
};

// Cox-Ross-Rubinstein with drift (CRR-drift)
class CRR_drift: public BinomialStrategy
{
public:
    double eta;
public:
    CRR_drift(double sigma, double r, double dt, double _eta);
    virtual ~CRR_drift();
};

// Leisen-Reimer (LR)
class LeisenReimer: public BinomialStrategy
{
public:
    double p_bar;
    int n;
    double S0;
    double K;
    double T;
public:
    LeisenReimer(double sigma, double r, double dt, double _S0, double _K, double _T, int _n);
    virtual ~LeisenReimer();
};

#endif // BINOMIALSTRATEGY_H
