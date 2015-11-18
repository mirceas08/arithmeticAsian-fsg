#include <iostream>
#include <math.h>
#include <algorithm>
#include <string>
#include "forward_shooting_grid.h"

#include <armadillo>
using namespace arma;

int main(int argc, char** argv)
{
    std::string dataFile = "data/data.dat";

    cout << "Option price: " << forwardShootingGrid(dataFile.c_str());
    return 0;
}
