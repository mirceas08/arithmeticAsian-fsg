#include <armadillo>
using namespace arma;

void interpolatePrices(const vec &av, const vec &optionPrice, const vec &F, vec &interpolatedValue)
{
    int vecSize = interpolatedValue.size();
    int spot;
    double diff;

    for (int k = 0; k < vecSize; k++) {
        uvec spotVec = find(av <= F(k), 1, "last");
        spot = as_scalar(spotVec);

        if (spotVec.size() == 0) {
            interpolatedValue(k) = optionPrice(0);
        }
        else if (spot == vecSize-1) {
            interpolatedValue(k) = optionPrice(vecSize-1);
        }
        else
        {
            diff = av(spot+1) - av(spot);

            interpolatedValue(k) = ( F(k) - av(spot) ) * optionPrice(spot+1) + ( av(spot+1) - F(k) ) * optionPrice(spot);
            interpolatedValue(k) /= diff;
        }
    }
}
