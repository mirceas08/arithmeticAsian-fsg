#include <math.h>
#include <vector>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>

#include "helpers.h"
#include "binomialStrategy.h"

#include <armadillo>
using namespace arma;

double fsg_vector(std::string dataFile)
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
    std::string treeStrategy;
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
        else if (variable == "treeStrategy")
            treeStrategy = value;
        else
            break;
    }

    transform(putCall.begin(), putCall.end(), putCall.begin(), ::toupper);
    transform(optionStyle.begin(), optionStyle.end(), optionStyle.begin(), ::toupper);
    double dt = T / static_cast<double>(n);

    BinomialStrategy* latticeStrategy;
    if (treeStrategy == "CRR")
        latticeStrategy = new CRR(sigma, r, dt);
    else if (treeStrategy == "JR")
        latticeStrategy = new JarrowRudd(sigma, r, dt);


    double u = latticeStrategy->u;
    double d = latticeStrategy->d;
    double p = (std::exp(latticeStrategy->r * latticeStrategy->dt) - latticeStrategy->d) / (latticeStrategy->u - latticeStrategy->d);
    double q = 1.0 - p;

//    double dt = T / static_cast<double>(n);
//    double u = std::exp(sigma * std::sqrt(dt));
//    double d = std::exp(-sigma * std::sqrt(dt));
//    double p = (std::exp(r*dt) - d) / (u-d);
//    double q = 1.0 - p;

    /* ------------------------ Initialize some values ------------------------ */
    vec strike(numAverages);
    strike.fill(K);
    vec zeroVec = zeros<vec>(numAverages);

    mat S(n+1, n+1);
    std::vector<field<vec> > av;
    std::vector<vec> optionPrice;

    /* ------------------------ Build binomial tree for S ------------------------ */
	for (j = 0; j <= n; j++) {
        for (i = 0; i <= j; i++) {
            S(i,j) = S0 * pow(u, j-i) * pow(d, i);
        }
	}

	/* ------------------------ Shoot averages ------------------------ */
	for (j = 0; j <= n; j++) {
        field<vec> temp_field(j+1);

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
            temp_field(i) = temp;
        }
        av.push_back(temp_field);
	}

	/* ------------------------ Compute terminal payoffs ------------------------ */
	vec temp_optionPrice(numAverages);
	for (i = 0; i <= n; i++) {
        temp_optionPrice = max(av[n].at(i) - strike, zeroVec);
        optionPrice.push_back(temp_optionPrice);
	}

	/* ------------------------ Backward recursion ------------------------ */
    for (j = n-1; j >= 0; j--) {
        std::vector<vec> optionPrice_temp;
        for (i = 0; i <= j; i++) {
            vec currentS(numAverages);
            currentS.fill(S(i,j));
            vec Au = ((j+1) * av[j].at(i) + u * currentS) / (j + 2);
            vec Ad = ((j+1) * av[j].at(i) + d * currentS) / (j + 2);

            vec interpOption_u(numAverages), interpOption_d(numAverages);
            interpolatePrices(av[j+1].at(i), optionPrice[i], Au, interpOption_u);
            interpolatePrices(av[j+1].at(i+1), optionPrice[i+1], Ad, interpOption_d);

            vec temp_optionPrice(numAverages);
            if (optionStyle == "E") {
                temp_optionPrice = std::exp(-r*dt) * ( p * interpOption_u + q * interpOption_d );
            }
            else if (optionStyle == "A") {
                vec continuationValue = std::exp(-r*dt) * ( p * interpOption_u + q * interpOption_d );
                vec exerciseValue = max(av[j].at(i) - strike, zeroVec);
                temp_optionPrice = max(continuationValue, exerciseValue);
            }
            else {
                cout << "optionStyle can be either E (european) or A (american)" << endl;
                return -1;
            }
            optionPrice_temp.push_back(temp_optionPrice);
        }

        optionPrice = optionPrice_temp;
    }

    delete latticeStrategy;
    //return option price
    return optionPrice[0].at(0);
}

