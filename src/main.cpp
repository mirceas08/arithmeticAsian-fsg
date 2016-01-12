#ifndef ARITHMETICAVERAGE_FSG_CPP
#define ARITHMETICAVERAGE_FSG_CPP

#include <iostream>
#include <math.h>
#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>
#include "fsg_average.h"
#include "binomialStrategy.h"
#include "option.h"

#include <armadillo>
using namespace arma;

int main(int argc, char** argv)
{
    wall_clock timer;
    timer.tic();

    std::string dataFile = argv[1]; // get data file from command line
    cout << "Option price: " << fsg_average(dataFile) << std::endl; // output option price

    // timer
    double timeElapsed = timer.toc();
    cout << "Time elapsed: " << timeElapsed << std::endl;

    return 0;
}

#endif // ARITHMETICAVERAGE_FSG_CPP
