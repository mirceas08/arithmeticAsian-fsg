#include <iostream>
#include <math.h>
#include <algorithm>

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

    double dt = T / static_cast<double>(n);
    double u = std::exp(sigma * std::sqrt(dt));
    double d = 1.0 / u;
    double p = (std::exp(r*dt) - d) / (u-d);

    mat S(n+1, n+1);
    field<vec> av(n+1, n+1);
    av.fill(zeros<vec>(numAverages));

    // Build the binomial tree for the stock
	for (j = 0; j <= n; j++) {
        for (i = 0; i <= j; i++) {
            S(i,j) = S0 * pow(u, j-i) * pow(d, i);
        }
	}

    // Build the binomial tree for averages
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

    cout << S(2,4) << endl;
    cout << av(2,4) << endl;

    return 0;
}
