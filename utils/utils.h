//
// Created by Thomas Gueutal on 22/02/2021.
//

#ifndef ENGINE_UTILS_H
#define ENGINE_UTILS_H


/*
 * Round a double to an int, using std::round(double)
 *
 * @param d:
 *      The value to be rounded.
 */
int roundToInt(double d);


/*
 * Convert degrees to radians
 *
 * @param degrees:
 *      The value in degrees to convert to radian
 */
double toRadians(double degrees);


/*
 * Clamp the passed value according to the given bounds.
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
unsigned int clamp(unsigned int val, unsigned int lowerBound, unsigned int upperBound);


#endif //ENGINE_UTILS_H
