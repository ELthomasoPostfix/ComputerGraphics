//
// Created by Thomas Gueutal on 13/02/2021.
//

#include "utils.h"
#include <cmath>



int roundToInt(double d)
{
    return static_cast<int>(std::round(d));
}


double toRadians(const double degrees) {
    return ( degrees * M_PI ) / 180.0;
}


unsigned int clamp(unsigned int val, unsigned int lowerBound, unsigned int upperBound) {
    if (val < lowerBound)
        return lowerBound;
    else if (val > upperBound)
        return upperBound;

    return val;
}
