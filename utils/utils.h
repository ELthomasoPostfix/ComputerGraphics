//
// Created by Thomas Gueutal on 22/02/2021.
//

#ifndef ENGINE_UTILS_H
#define ENGINE_UTILS_H

#include <string>

/*
 * Description:
 *      Round a double to an int, using std::round(double)
 *
 * @param d:
 *      The value to be rounded.
 */
int roundToInt(double d);


/*
 * Description:
 *      Convert degrees to radians
 *
 * @param degrees:
 *      The value in degrees to convert to radian
 */
double toRadians(double degrees);


/*
 * Description:
 *      Clamp the passed value according to the given bounds.
 *      * val < lowerBound returns lowerBound
 *      * val > upperBound returns upperBound
 *      * else returns val
 *
 * @param val:
 *      The value to be clamped.
 * @param lowerBound:
 *      If val < lowerBound, clamp val to lowerBound.
 * @param upperBound:
 *      If val > upperBound, clamp val to upperBound.
 */
int clamp(int val, int lowerBound, int upperBound);
double clamp(double val, double lowerBound, double upperBound);

/*
 * Description:
 *      We assume the filePath argument is the path of to a file.
 *      We truncate the file name (including extension) until either
 *      the first '/' from the end encountered or until the result
 *      is an empty string.
 *
 * @param filePath:
 *      The file path to truncate.
 */
std::string truncateFileName(const std::string& filePath);



#endif //ENGINE_UTILS_H
