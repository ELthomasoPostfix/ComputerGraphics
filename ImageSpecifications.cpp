//
// Created by Thomas Gueutal on 10/05/2021.
//

#include "ImageSpecifications.h"


ImageSpecifications::ImageSpecifications(const L2D::Lines2D &lines, const double size) {


    double minX =  std::numeric_limits<double>::infinity();
    double minY =  std::numeric_limits<double>::infinity();
    double maxX = -std::numeric_limits<double>::infinity();
    double maxY = -std::numeric_limits<double>::infinity();

    determineExtrema(lines, minX, maxX, minY, maxY);

    // find the differences between the smallest and largest x and y
    double rangeX = maxX - minX;

    double rangeY = maxY - minY;


    // The image should have at least a width or a height different from 0
    // A width AND height of 0 would imply no image can be generated
    assert(std::abs(rangeX) > 0.0 || std::abs(rangeY) > 0.0);

    // if rangeX or rangeY is 0, then set it to 1
    // ==> ensure that the width and height of the image will both be > 0
    rangeX = (rangeX > 0.0) ? rangeX : 1.0;
    rangeY = (rangeY > 0.0) ? rangeY : 1.0;


    // one of the images will be assigned size and the other a percent (range/max(ranges)) of size
    this->imageX = size * (rangeX / std::max(rangeX, rangeY));

    this->imageY = size * (rangeY / std::max(rangeX, rangeY));

    // find the factor by which to scale every line
    this->projectionScreenDistance = 0.95 * (this->imageX / rangeX);

    // Calculate the rescaled center of the set of lines
    double DCx = this->projectionScreenDistance * ((minX + maxX) / 2.0);
    double DCy = this->projectionScreenDistance * ((minY + maxY) / 2.0);

    // Find the new center of the lines
    this->dx = (this->imageX / 2.0) - DCx;
    this->dy = (this->imageY / 2.0) - DCy;
}

ImageSpecifications::ImageSpecifications(const L2D::Lines2DZ &lines, const double size) {


    double minX =  std::numeric_limits<double>::infinity();
    double minY =  std::numeric_limits<double>::infinity();
    double maxX = -std::numeric_limits<double>::infinity();
    double maxY = -std::numeric_limits<double>::infinity();

    determineExtrema(lines, minX, maxX, minY, maxY);

    // find the differences between the smallest and largest x and y
    double rangeX = maxX - minX;

    double rangeY = maxY - minY;


    // The image should have at least a width or a height different from 0
    // A width AND height of 0 would imply no image can be generated
    assert(std::abs(rangeX) > 0.0 || std::abs(rangeY) > 0.0);

    // if rangeX or rangeY is 0, then set it to 1
    // ==> ensure that the width and height of the image will both be > 0
    rangeX = (rangeX > 0.0) ? rangeX : 1.0;
    rangeY = (rangeY > 0.0) ? rangeY : 1.0;


    // one of the images will be assigned size and the other a percent (range/max(ranges)) of size
    this->imageX = size * (rangeX / std::max(rangeX, rangeY));

    this->imageY = size * (rangeY / std::max(rangeX, rangeY));

    // find the factor by which to scale every line
    this->projectionScreenDistance = 0.95 * (this->imageX / rangeX);

    // Calculate the rescaled center of the set of lines
    double DCx = this->projectionScreenDistance * ((minX + maxX) / 2.0);
    double DCy = this->projectionScreenDistance * ((minY + maxY) / 2.0);

    // Find the new center of the lines
    this->dx = (this->imageX / 2.0) - DCx;
    this->dy = (this->imageY / 2.0) - DCy;
}

ImageSpecifications::ImageSpecifications(const L3D::Figures3D &figures, const double size) {

    double minX =  std::numeric_limits<double>::infinity();
    double minY =  std::numeric_limits<double>::infinity();
    double maxX = -std::numeric_limits<double>::infinity();
    double maxY = -std::numeric_limits<double>::infinity();

    determineExtrema(figures, minX, maxX, minY, maxY);

    // find the differences between the smallest and largest x and y
    double rangeX = maxX - minX;

    double rangeY = maxY - minY;


    // The image should have at least a width or a height different from 0
    // A width AND height of 0 would imply no image can be generated
    assert(std::abs(rangeX) > 0.0 || std::abs(rangeY) > 0.0);

    // if rangeX or rangeY is 0, then set it to 1
    // ==> ensure that the width and height of the image will both be > 0
    rangeX = (rangeX > 0.0) ? rangeX : 1.0;
    rangeY = (rangeY > 0.0) ? rangeY : 1.0;


    // one of the images will be assigned size and the other a percent (range/max(ranges)) of size
    this->imageX = size * (rangeX / std::max(rangeX, rangeY));

    this->imageY = size * (rangeY / std::max(rangeX, rangeY));

    // find the factor by which to scale every line
    this->projectionScreenDistance = 0.95 * (this->imageX / rangeX);

    // Calculate the rescaled center of the set of lines
    double DCx = this->projectionScreenDistance * ((minX + maxX) / 2.0);
    double DCy = this->projectionScreenDistance * ((minY + maxY) / 2.0);

    // Find the new center of the lines
    this->dx = (this->imageX / 2.0) - DCx;
    this->dy = (this->imageY / 2.0) - DCy;
}

void ImageSpecifications::determineExtrema(const L2D::Lines2D &lines, double &minX, double &maxX, double &minY, double &maxY) {

    // variables to classify the coordinates of any given line
    double smallerX, largerX, smallerY, largerY;

    for (const L2D::Line2D& line : lines) {
        smallerX = std::min(line.p1.x, line.p2.x);
        largerX = std::max(line.p1.x, line.p2.x);
        smallerY = std::min(line.p1.y, line.p2.y);
        largerY = std::max(line.p1.y, line.p2.y);

        if (smallerX < minX)
            minX = smallerX;
        if (largerX > maxX)
            maxX = largerX;
        if (smallerY < minY)
            minY = smallerY;
        if (largerY > maxY)
            maxY = largerY;
    }
}

void ImageSpecifications::determineExtrema(const L2D::Lines2DZ &lines, double &minX, double &maxX, double &minY, double &maxY) {

    // variables to classify the coordinates of any given line
    double smallerX, largerX, smallerY, largerY;

    for (const L2D::Line2DZ& line : lines) {
        smallerX = std::min(line.p1.x, line.p2.x);
        largerX = std::max(line.p1.x, line.p2.x);
        smallerY = std::min(line.p1.y, line.p2.y);
        largerY = std::max(line.p1.y, line.p2.y);

        if (smallerX < minX)
            minX = smallerX;
        if (largerX > maxX)
            maxX = largerX;
        if (smallerY < minY)
            minY = smallerY;
        if (largerY > maxY)
            maxY = largerY;
    }
}

void ImageSpecifications::determineExtrema(const L3D::Figures3D &figures,
                                           double &minX, double &maxX,
                                           double &minY, double &maxY) {

    // variables to classify the coordinates of any given line
    L2D::Point2D projectedP;

    for (const L3D::Figure& fig : figures) {
        for (const Vector3D& point3D : fig.points) {
            projectedP = L3D::projectPoint3D(point3D, 1.0);

            if (projectedP.x < minX)
                minX = projectedP.x;
            if (projectedP.x > maxX)
                maxX = projectedP.x;
            if (projectedP.y < minY)
                minY = projectedP.y;
            if (projectedP.y > maxY)
                maxY = projectedP.y;
        }
    }
}

L2D::Point2D ImageSpecifications::getMoveVector() const {
    return {this->dx, this->dy};
}
