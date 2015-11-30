#ifndef ARITHMETICAVERAGE_FSG_CPP
#define ARITHMETICAVERAGE_FSG_CPP

#include <iostream>
#include <math.h>
#include <algorithm>
#include <string>
#include "forwardShootingGrid.h"
#include "binomialStrategy.h"
#include "option.h"

#include <armadillo>
using namespace arma;

int main(int argc, char** argv)
{
    wall_clock timer;
    timer.tic();

    std::string dataFile = "data/data.dat";

    cout << "Option price: " << forwardShootingGrid(dataFile) << std::endl;

    double timeElapsed = timer.toc();
    cout << "Time elapsed: " << timeElapsed << std::endl;

    return 0;
}

#endif // ARITHMETICAVERAGE_FSG_CPP
