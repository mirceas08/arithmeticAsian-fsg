#ifndef HELPERS_H
#define HELPERS_H

#include <armadillo>
using namespace arma;

// function for interpolation of vectors that handles out-of-range values
// the built-in function interp1() of Armadillo doesn't deal with values out of range
void interpolate(const vec &X, const vec &Y, const vec &XX, vec &YY, std::string interpolationType)
{
    int vecSize = YY.size();
    int spot;

    if (interpolationType == "linear") {
        for (int k = 0; k < vecSize; k++) {
            uvec spotVec = find(X <= XX(k), 1, "last"); // find lower index
            spot = as_scalar(spotVec);

            if (spotVec.size() == 0) {
                YY(k) = Y(0);  // if value to interpolate is lower than the min of av take the ceiling
            }
            else if (spot == vecSize-1) {
                YY(k) = Y(vecSize-1); // if value to interpolate is greater than the max of av take the floor
            }
            else { // otherwise do linear interpolation
                YY(k) = Y(spot) + ( Y(spot+1) - Y(spot) ) * ( XX(k) - X(spot) ) / ( X(spot+1) - X(spot) );
            }
        }
    }
    // performs badly
    else if (interpolationType == "floor") {
        for (int k = 0; k < vecSize; k++) {
            uvec spotVec = find(X <= XX(k), 1, "last"); // find lower index
            spot = as_scalar(spotVec);

            if (spotVec.size() == 0) {
                YY(k) = Y(0);  // if value to interpolate is lower than the min of av take the ceiling
            }
            else { // otherwise do linear interpolation
                YY(k) = Y(spot);
            }
        }
    }
    // performs badly
    else if (interpolationType == "ceil") {
        for (int k = 0; k < vecSize; k++) {
            uvec spotVec = find(X <= XX(k), 1, "last"); // find lower index
            spot = as_scalar(spotVec);

            if (spotVec.size() == 0) {
                YY(k) = Y(0);  // if value to interpolate is lower than the min of av take the ceiling
            }
            else if (spot == vecSize-1) {
                YY(k) = Y(vecSize-1); // if value to interpolate is greater than the max of av take the floor
            }
            else { // otherwise do linear interpolation
                YY(k) = Y(spot+1);
            }
        }
    }
    /*
    the version for quadratic interpolation is commented out because it is buggy along the edges of the tree
    linear interpolation performs well anyway
    */
//    else if (interpolationType == "quadratic") {
//        for (int k = 0; k < vecSize; k++) {
//            uvec spotVec = find(X <= XX(k), 1, "last"); // find lower index
//            spot = as_scalar(spotVec);
//            double b0, b1, b2;
//
//            if (spotVec.size() == 0) {
//                YY(k) = Y(0);
//            }
//            else if (spot == vecSize-1) {
//                YY(k) = Y(vecSize-1); // if value to interpolate is greater than the max of av take the floor
//            }
//            else if (spot == vecSize-2 && spot == vecSize-3) {
//                b0 = Y(spot-1);
//                b1 = ( Y(spot) - Y(spot-1) ) / ( X(spot) - X(spot-1) );
//                b2 = ( ( Y(spot+1) - Y(spot) ) / ( X(spot+1) - X(spot) ) - ( Y(spot) - Y(spot-1) ) / ( X(spot) - X(spot-1) ) ) / ( X(spot+1) - X(spot-1) );
//
//                YY(k) = b0 + b1 * ( XX(k) - X(spot-1) ) + b2 * ( XX(k) - X(spot-1) ) * ( XX(k) - X(spot) );
//            }
//            else { // otherwise do quadratic interpolation
//                b0 = Y(spot);
//                b1 = ( Y(spot+1) - Y(spot) ) / ( X(spot+1) - X(spot) );
//                b2 = ( ( Y(spot+2) - Y(spot+1) ) / ( X(spot+2) - X(spot+1) ) - ( Y(spot+1) - Y(spot) ) / ( X(spot+1) - X(spot) ) ) / ( X(spot+2) - X(spot) );
//
//                YY(k) = b0 + b1 * ( XX(k) - X(spot) ) + b2 * ( XX(k) - X(spot) ) * ( XX(k) - X(spot+1) );
//            }
//        }
//    }
    else { // linear interpolation as default
        for (int k = 0; k < vecSize; k++) {
            uvec spotVec = find(X <= XX(k), 1, "last"); // find lower index
            spot = as_scalar(spotVec);

            if (spotVec.size() == 0) {
                YY(k) = Y(0);  // if value to interpolate is lower than the min of av take the ceiling
            }
            else if (spot == vecSize-1) {
                YY(k) = Y(vecSize-1); // if value to interpolate is greater than the max of av take the floor
            }
            else {  // otherwise do linear interpolation
                YY(k) = Y(spot) + ( Y(spot+1) - Y(spot) ) * ( XX(k) - X(spot) ) / ( X(spot+1) - X(spot) );
            }
        }
    }
}

// function to generate a vector with values (linearly or logarithmically spaced) within a range
// logspace built-in function not available in Armadillo
vec spacedVector(const double &a, const double &b, int N, std::string spaceType)
{
    vec myVector(N);

    if (spaceType == "linspace") {
        double diff = b - a;
        double spacing = diff / (N-1);

        for (int k = 0; k < N; k++) {
            myVector(k) = a + k * spacing;
        }
    }
    else if (spaceType == "logspace") {
        double logdiff = std::log(b) - std::log(a);
        double logspacing = logdiff / (N-1);

        for (int k = 0; k < N; k++) {
            myVector(k) = std::exp(std::log(a) + k * logspacing);
        }
    }
    else {
        double diff = b - a;
        double spacing = diff / (N-1);

        for (int k = 0; k < N; k++) {
            myVector(k) = a + k * spacing;
        }
    }

    return myVector;
}

// friend function for Leisen-Reimer class
double h(double x, int n)
{
    double sign;
    if (x > 0)
        sign = 1.0;
    else if (x < 0)
        sign = -1.0;
    else
        sign = 0.0;

    double argument;
    double parenthesis;

    double term = (x / (n + 1/3 + 0.1 / (n+1)));
    parenthesis = -1*(n + 1/6) * term * term;
    return 0.5 + sign * std::sqrt(0.25 - 0.25 * std::exp(parenthesis));
}


// CDF of standard normal distribution
double normalCDF(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}

#endif // HELPERS_H
