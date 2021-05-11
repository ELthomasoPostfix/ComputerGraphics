//
// Created by Thomas Gueutal on 10/05/2021.
//

#ifndef ENGINE_IMAGESPECIFICATIONS_H
#define ENGINE_IMAGESPECIFICATIONS_H

#include "L3D.h"


/*
 * Description:
 *      This struct computes computes and stores the specifications of the img::EasyImage
 *      that is to be generated. This is done based on either a list of projected 2D lines
 *      (L2D::Line2D or L2D::Line2DZ) or a L3D::Figures3D object.
 */
struct ImageSpecifications {
    public:
        /*
         * Description:
         *      A constructor meant to be used for extracting the img::EasyImage
         *      specifications from a list of L2D::Lines2D or L2D::Lines2DZ.
         *
         * @param lines:
         *      The list of lines from which to extract the img::EasyImage
         *      specifications.
         * @param size:
         *      The size that the longest side of the image should be.
         */
        ImageSpecifications(const L2D::Lines2D& lines, double size);
        ImageSpecifications(const L2D::Lines2DZ& lines, double size);

        /*
         * Description:
         *      A constructor meant to be used for extracting the img::EasyImage
         *      specifications from a L3D::Figures3D list.
         *
         * @param figures:
         *      The list of L3D::Figure's from which to extract the img::EasyImage
         *      specifications.
         * @param size:
         *      The size that the longest side of the image should be.
         */
        ImageSpecifications(const L3D::Figures3D& figures, double size);

        /*
         * Description:
         *      Determine which vector should be added to each of the projected points
         *      in order to fit all projected lines onto the to generate image.
         */
        L2D::Point2D getMoveVector() const;

    private:
        static void determineExtrema(const L2D::Lines2D& lines,
                              double& minX, double& maxX,
                              double& minY, double& maxY);

        static void determineExtrema(const L2D::Lines2DZ& lines,
                              double& minX, double& maxX,
                              double& minY, double& maxY);

        static void determineExtrema(const L3D::Figures3D & figures,
                              double& minX, double& maxX,
                              double& minY, double& maxY);

    public:
        double imageX;
        double imageY;
        double dx;
        double dy;
        double projectionScreenDistance;
};



#endif //ENGINE_IMAGESPECIFICATIONS_H
