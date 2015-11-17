#include <armadillo>
using namespace arma;

void interpolatePrices(const vec &av, const vec &F, const vec &optionPrice, vec &interpolatedValue)
{
    int vecSize = interpolatedValue.size();
    int spot;
    double diff;

    for (int k = 0; k < vecSize; k++) {
        spot = as_scalar(find(av <= F(k), 1, "last"));
        diff = av(spot+1) - av(spot);

        interpolatedValue(k) = ( F(k) - av(spot) ) * optionPrice(spot+1) + ( av(spot+1) - F(k) ) * optionPrice(spot);
        interpolatedValue(k) /= diff;
    }
}
