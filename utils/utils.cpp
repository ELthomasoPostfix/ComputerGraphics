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


int clamp(const int val, const int lowerBound, const int upperBound) {
    if (val < lowerBound)
        return lowerBound;
    else if (val > upperBound)
        return upperBound;

    return val;
}

double clamp(const double val, const double lowerBound, const double upperBound) {
    if (val < lowerBound)
        return lowerBound;
    else if (val > upperBound)
        return upperBound;

    return val;
}


std::string truncateFileName(const std::string& filePath) {
    std::string truncatedPath = filePath;

    // filePath.size() returns and unsigned type, so size - 1 could wrap around --> cast to int
    int i = ((int) filePath.size()) - 1;

    while (i >= 0) {

        if (filePath.at(i) == '/')
            break;

        i--;
    }

    truncatedPath.resize(i+1);

    return truncatedPath;
}
