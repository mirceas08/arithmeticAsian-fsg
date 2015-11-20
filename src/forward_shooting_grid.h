#ifndef FORWARDSHOOTINGGRID_H
#define FORWARDSHOOTINGGRID_H

#include <math.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>
#include "helpers.h"

#include <armadillo>
using namespace arma;

double forwardShootingGrid(std::string dataFile)
{
    int i, j;

    /* ------------------------ Set starting values from file ------------------------ */
    int n;
    double S0;
    double K;
    double r;
    double T;
    double sigma;
    int numAverages;
    std::string putCall;
    std::string optionStyle;

    std::ifstream fIN(dataFile.c_str());
    std::string line;
    getline(fIN, line); // skip header

    while (std::getline(fIN, line)) {
        std::stringstream stream(line);
        std::string variable;
        std::string value;

        stream >> variable >> value;

        if (variable == "nsteps")
            n = atof(value.c_str());
        else if (variable == "numAverages")
            numAverages = atof(value.c_str());
        else if (variable == "spot")
            S0 = atof(value.c_str());
        else if (variable == "strike")
            K = atof(value.c_str());
        else if (variable == "r")
            r = atof(value.c_str());
        else if (variable == "maturity")
            T = atof(value.c_str());
        else if (variable == "sigma")
            sigma = atof(value.c_str());
        else if (variable == "putCall")
            putCall = value;
        else if (variable == "optionStyle")
            optionStyle = value;
        else
            break;
    }

    transform(putCall.begin(), putCall.end(), putCall.begin(), ::toupper);
    transform(optionStyle.begin(), optionStyle.end(), optionStyle.begin(), ::toupper);

    /* ------------------------ Initialize some values ------------------------ */
    vec strike(numAverages);
    strike.fill(K);
    vec zeroVec = zeros<vec>(numAverages);

    double dt = T / static_cast<double>(n);
    double u = std::exp(sigma * std::sqrt(dt));
    double d = std::exp(-sigma * std::sqrt(dt));
    double p = (std::exp(r*dt) - d) / (u-d);
    double q = 1.0 - p;

    sp_mat S(n+1, n+1);
    field<vec> av(n+1, n+1);
    field<vec> optionPrice(n+1);

    /* ------------------------ Build binomial tree for S ------------------------ */
	for (j = 0; j <= n; j++) {
        for (i = 0; i <= j; i++) {
            S(i,j) = S0 * pow(u, j-i) * pow(d, i);
        }
	}

	/* ------------------------ Shoot averages ------------------------ */
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

            vec temp(numAverages);
            for (int k = 0; k < numAverages; k++) {
                temp(k) = a_min + k * spacing;
            }
            av(i,j) = temp;
        }
	}

	/* ------------------------ Compute terminal payoffs ------------------------ */
	for (i = 0; i <= n; i++) {
        optionPrice(i) = max(av(i,n) - strike, zeroVec);
	}

	/* ------------------------ Backward recursion ------------------------ */
    for (j = n-1; j >= 0; j--) {
        for (i = 0; i <= j; i++) {
            vec currentS(numAverages);
            currentS.fill(S(i,j));
            vec Au = ((j+1) * av(i,j) + u * currentS) / (j + 2);
            vec Ad = ((j+1) * av(i,j) + d * currentS) / (j + 2);

            vec interpOption_u(numAverages), interpOption_d(numAverages);
            interpolatePrices(av(i,j+1), optionPrice(i), Au, interpOption_u);
            interpolatePrices(av(i+1,j+1), optionPrice(i+1), Ad, interpOption_d);

            if (optionStyle == "E") {
                optionPrice(i) = std::exp(-r*dt) * ( p * interpOption_u + q * interpOption_d );
            }
            else if (optionStyle == "A") {
                vec continuationValue = std::exp(-r*dt) * ( p * interpOption_u + q * interpOption_d );
                vec exerciseValue = max(av(i,j) - strike, zeroVec);
                optionPrice(i) = max(continuationValue, exerciseValue);
            }
            else {
                cout << "optionStyle can be either E (european) or A (american)" << endl;
                return -1;
            }

        }
    }

    //return option price
    return optionPrice(0).at(0);
}

#endif // FORWARDSHOOTINGGRID_H
