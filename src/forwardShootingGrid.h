#include <math.h>
#include <vector>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>

#include "helpers.h"
#include "binomialStrategy.h"
#include "option.h"

#include <armadillo>
using namespace arma;

double forwardShootingGrid(std::string dataFile)
{
    int i, j; // counters

    /* ------------------------ Set starting values from file ------------------------ */
    int n;
    double S0;
    double K;
    double r;
    double T;
    double sigma;
    int numAverages;
    std::string putCall;
    std::string optionType;
    std::string optionStyle;
    std::string treeStrategy;
    std::string interpolationType;
    std::string spaceType;

    std::ifstream fIN(dataFile.c_str());
    std::string line;

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
        else if (variable == "optionType")
            optionType = value;
        else if (variable == "optionStyle")
            optionStyle = value;
        else if (variable == "treeStrategy")
            treeStrategy = value;
        else if (variable == "interpolationType")
            interpolationType = value;
        else if (variable == "spaceType")
            spaceType = value;
    }

    transform(putCall.begin(), putCall.end(), putCall.begin(), ::toupper);
    transform(optionType.begin(), optionType.end(), optionType.begin(), ::toupper);
    transform(optionStyle.begin(), optionStyle.end(), optionStyle.begin(), ::toupper);
    double dt = T / static_cast<double>(n);

    /* ------------------------ Set tree parameters according to the lattice strategy ------------------------ */

    BinomialStrategy* latticeStrategy; // latticeStrategy allocated on the heap
    if (treeStrategy == "CRR")
        latticeStrategy = new CRR(sigma, r, dt);
    else if (treeStrategy == "JR")
        latticeStrategy = new JarrowRudd(sigma, r, dt);
    else
        latticeStrategy = new CRR(sigma, r, dt); // default strategy


    double u = latticeStrategy->u;
    double d = latticeStrategy->d;
    double p = (std::exp(latticeStrategy->r * latticeStrategy->dt) - d) / (u-d);
    double q = 1.0 - p;

    /* ------------------------ Initialize some values ------------------------ */
    vec strike(numAverages);
    strike.fill(K);
    vec zeroVec = zeros<vec>(numAverages);

    sp_mat S(n+1, n+1);
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

            vec temp(numAverages);
            spacedVector(a_min, a_max, temp, spaceType);
            temp_field(i) = temp;
        }
        av.push_back(temp_field);
	}

    /* ------------------------ Construct option object ------------------------ */
    Option* pathOption;
    if (putCall == "CALL")
        pathOption = new AverageCallOption(strike);
    else if (putCall == "PUT")
        pathOption = new AveragePutOption(strike);
    else {
        cout << "!!! putCall parameter can be either call or put" << endl;
        return -1;
    }

	/* ------------------------ Compute terminal payoffs ------------------------ */
	vec temp_optionPrice(numAverages);
//	for (i = 0; i <= n; i++) {
//        temp_optionPrice = max(av[n].at(i) - strike, zeroVec);
//        optionPrice.push_back(temp_optionPrice);
//	}

   for (i = 0; i <= n; i++) {
        temp_optionPrice = pathOption->payoff(av[n].at(i));
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
            interpolate(av[j+1].at(i), optionPrice[i], Au, interpOption_u, interpolationType);
            interpolate(av[j+1].at(i+1), optionPrice[i+1], Ad, interpOption_d, interpolationType);

            vec temp_optionPrice(numAverages);
            if (optionStyle == "E") {
                temp_optionPrice = std::exp(-r*dt) * ( p * interpOption_u + q * interpOption_d );
            }
            else if (optionStyle == "A") {
                vec continuationValue = std::exp(-r*dt) * ( p * interpOption_u + q * interpOption_d );
                vec exerciseValue = pathOption->payoff(av[j].at(i));
                temp_optionPrice = max(continuationValue, exerciseValue);
            }
            else {
                cout << "!!! optionStyle parameter can be either E (european) or A (american)" << endl;
                return -1;
            }
            optionPrice_temp.push_back(temp_optionPrice);
        }

        optionPrice = optionPrice_temp;
    }

    /* ------------------------ Clean from the heap and return ------------------------ */

    delete latticeStrategy;
    delete pathOption;
    return optionPrice[0].at(0);
}

