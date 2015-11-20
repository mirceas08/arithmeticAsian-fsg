#ifndef ARITHMETICAVERAGE_FSG_CPP
#define ARITHMETICAVERAGE_FSG_CPP

#include <iostream>
#include <math.h>
#include <algorithm>
#include <string>
#include "forward_shooting_grid.h"
#include "fsg_vector.h"

#include <armadillo>
using namespace arma;

int main(int argc, char** argv)
{
    std::string dataFile = "data/data.dat";

    cout << "Option price: " << forwardShootingGrid(dataFile.c_str()) << std::endl;
    //cout << "Option price with <vector>: " << fsg_vector(dataFile.c_str()) << std::endl;
    return 0;
}

#endif // ARITHMETICAVERAGE_FSG_CPP
