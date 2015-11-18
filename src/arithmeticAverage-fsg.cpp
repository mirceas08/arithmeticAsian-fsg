#include <iostream>
#include <math.h>
#include <algorithm>
#include <string>
#include "forward_shooting_grid.h"

#include <armadillo>
using namespace arma;

int main(int argc, char** argv)
{
    int i, j;

    std::string optionStyle = "A";
    int n = 60;
    int numAverages = 100;
    double S0 = 50.0;
    double sigma = 0.40;
    double r = 0.10;
    double T = 1.0;
    double K = 50.0;

    cout << "Option price: " << forwardShootingGrid(n, S0, K, r, sigma, T, numAverages, optionStyle);
    return 0;
}
