#ifndef HELPERS_H
#define HELPERS_H

#include <armadillo>
using namespace arma;

// function for interpolation of vectors that handles out-of-range values
// the built-in function interp1() of Armadillo doesn't deal with values out of range
void interpolate(const vec &X, const vec &Y, const vec &XX, vec &YY, std::string interpolationType)
{
    int vecSize = YY.size();

    if (interpolationType == "linear") {
        int spot;
        double diff;

        for (int k = 0; k < vecSize; k++) {
            uvec spotVec = find(X <= XX(k), 1, "last"); // find lower index
            spot = as_scalar(spotVec);

            if (spotVec.size() == 0) {
                YY(k) = Y(0);  // if value to interpolate is lower than the min of av take the ceiling
            }
            else if (spot == vecSize-1) {
                YY(k) = Y(vecSize-1); // if value to interpolate is greater than the max of av take the floor
            }
            else // otherwise do linear interpolation
            {
                diff = X(spot+1) - X(spot);

                YY(k) = ( XX(k) - X(spot) ) * Y(spot+1) + ( X(spot+1) - XX(k) ) * Y(spot);
                YY(k) /= diff;
            }
        }
    }
    else {
        int spot;
        double diff;

        for (int k = 0; k < vecSize; k++) {
            uvec spotVec = find(X <= XX(k), 1, "last"); // find lower index
            spot = as_scalar(spotVec);

            if (spotVec.size() == 0) {
                YY(k) = Y(0);  // if value to interpolate is lower than the min of av take the ceiling
            }
            else if (spot == vecSize-1) {
                YY(k) = Y(vecSize-1); // if value to interpolate is greater than the max of av take the floor
            }
            else // otherwise do linear interpolation
            {
                diff = X(spot+1) - X(spot);

                YY(k) = ( XX(k) - X(spot) ) * Y(spot+1) + ( X(spot+1) - XX(k) ) * Y(spot);
                YY(k) /= diff;
            }
        }
    }
}

// function to generate a vector with values (linearly or logarithmically spaced) within a range
// logspace built-in function not available in Armadillo
void spacedVector(const double &a, const double &b, vec &myVector, std::string spaceType)
{
    int vec_size = myVector.size();

    if (spaceType == "linspace") {
        double diff = b - a;
        double spacing = diff / (vec_size-1);

        for (int k = 0; k < vec_size; k++) {
            myVector(k) = a + k * spacing;
        }
    }
    else if (spaceType == "logspace") {
        double logdiff = std::log(b) - std::log(a);
        double logspacing = logdiff / (vec_size-1);

        for (int k = 0; k < vec_size; k++) {
            myVector(k) = std::exp(std::log(a) + k * logspacing);
        }
    }
    else {
        double diff = b - a;
        double spacing = diff / (vec_size-1);

        for (int k = 0; k < vec_size; k++) {
            myVector(k) = a + k * spacing;
        }
    }
}


#endif // HELPERS_H
