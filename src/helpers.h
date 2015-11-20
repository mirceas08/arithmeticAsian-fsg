#ifndef HELPERS_H
#define HELPERS_H

#include <armadillo>
using namespace arma;

// function for interpolation of values
void interpolatePrices(const vec &av, const vec &optionPrice, const vec &F, vec &interpolatedValue)
{
    int vecSize = interpolatedValue.size();
    int spot;
    double diff;

    for (int k = 0; k < vecSize; k++) {
        uvec spotVec = find(av <= F(k), 1, "last"); // find lower index
        spot = as_scalar(spotVec);

        if (spotVec.size() == 0) {
            interpolatedValue(k) = optionPrice(0);  // if value to interpolate is lower than the min of av take the ceiling
        }
        else if (spot == vecSize-1) {
            interpolatedValue(k) = optionPrice(vecSize-1); // if value to interpolate is greater than the max of av take the floor
        }
        else // otherwise do linear interpolation
        {
            diff = av(spot+1) - av(spot);

            interpolatedValue(k) = ( F(k) - av(spot) ) * optionPrice(spot+1) + ( av(spot+1) - F(k) ) * optionPrice(spot);
            interpolatedValue(k) /= diff;
        }
    }
}


#endif // HELPERS_H
