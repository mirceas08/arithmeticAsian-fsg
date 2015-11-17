#include <iostream>
#include <math.h>
#include <algorithm>
#include "helpers.h"

#include <armadillo>
using namespace arma;

int main(int argc, char** argv)
{
    int i, j;

    int n = 20;
    int numAverages = 4;
    double S0 = 50.0;
    double sigma = 0.40;
    double r = 0.10;
    double T = 1.0;
    double K = 50.0;
    vec strike(numAverages);
    strike.fill(K);
    vec zeroVec = zeros<vec>(numAverages);

    double dt = T / static_cast<double>(n);
    double u = std::exp(sigma * std::sqrt(dt));
    double d = std::exp(-sigma * std::sqrt(dt));
    double p = (std::exp(r*dt) - d) / (u-d);
    double q = 1.0 - p;

    mat S(n+1, n+1);
    field<vec> av(n+1, n+1);
    field<vec> interpolatedAV(n+1, n+1);
    field<vec> optionPrice(n+1, n+1);
    av.fill(zeros<vec>(numAverages));
    optionPrice.fill(zeros<vec>(numAverages));

    // Build the binomial tree for the stock
	for (j = 0; j <= n; j++) {
        for (i = 0; i <= j; i++) {
            S(i,j) = S0 * pow(u, j-i) * pow(d, i);
        }
	}

	// Shoot the averages
	for (j = 0; j <= n; j++) {
        for (i = 0; i <= j; i++) {
            double first_max = (1- pow(u, j - i + 1)) / (1-u);
            double second_max = pow(u, j-i) * d * (1 - pow(d, i)) / (1-d);

            double first_min = (1 - pow(d, i+1)) / (1-d);
            double second_min = pow(d, i) * u * (1 - pow(u, j - i)) / (1-u);

            double a_max = (S0 * first_max + S0 * second_max) / (j+1);
            double a_min = (S0 * first_min + S0 * second_min) / (j+1);
            double diff = a_max - a_min;
            double spacing = diff / (numAverages-1);

            for (int k = 0; k < numAverages; k++) {
                av(i,j)(k) = a_min + k * spacing;
            }
        }
	}

	// Compute terminal payoffs
	for (i = 0; i <= n; i++) {
        optionPrice(i,n) = max(av(i,n) - strike, zeroVec);
	}

    // Backward recursion
    for (j = n-1; j >= 0; j--) {
        for (i = 0; i <= j; i++) {
            vec currentS(numAverages);
            currentS.fill(S(i,j));
            vec Fu = ((j+1) * av(i,j) + u * currentS) / (j + 2);
            vec Fd = ((j+1) * av(i,j) + d * currentS) / (j + 2);

            vec interpOption_u(numAverages), interpOption_d(numAverages);
            interpolatePrices(av(i,j+1), optionPrice(i,j+1), Fu, interpOption_u);
            interpolatePrices(av(i+1,j+1), optionPrice(i+1,j+1), Fd, interpOption_d);

            optionPrice(i,j) = std::exp(-r*dt) * ( p * interpOption_u + q * interpOption_d );
        }
    }

    cout << "Option price: " << optionPrice(0,0)(0) << endl;

    return 0;
}
