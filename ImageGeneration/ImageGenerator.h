//
// Created by Thomas Gueutal on 10/05/2021.
//

#ifndef ENGINE_IMAGEGENERATOR_H
#define ENGINE_IMAGEGENERATOR_H

#include "ImageSpecifications.h"


namespace ImageGenerator {

    /*
     * Description:
     *      Given a ini::Configuration file, create and return all figures specified in
     *      the file in the form of a L3D::Figures3D object.
     *
     * @param configuration:
     *      The object representing the file which specifies all the needed
     *      L3D::Figure instances.
     * @param iniFilePath:
     *      The file path of the ini file. In the case that a lsystem should be
     *      generated, the L2D/L3D file should be located in the same directory as
     *      the ini file.
     * @param triangulate:
     *      Whether or not to triangulate any generated L3D::Figure.
     *      In the case that the "type" attribute of the "General" section is "ZBuffering",
     *      triangulate will be true.
     */
    L3D::Figures3D parseFigures(const ini::Configuration& configuration,
                                const std::string& iniFilePath = "",
                                bool triangulate = false);

    L3D::LightCaster parseLights(const ini::Configuration& configuration);

    /*
     * Description:
     *      project the L3D::Figures3D into a list of L2D::Lines2D.
     *      Does not make use ofz-buffering.
     *
     * @param figures:
     *      The figures to be projected.
     * @param projectionScreenDistance:
     *      The distance of the projections screen from the eye point.
     *      It is 1.0 by default.
     */
    L2D::Lines2D projectFigures(const L3D::Figures3D& figures, double projectionScreenDistance);

    /*
     * Description:
     *      project the L3D::Figures3D into a list of L2D::Lines2D.
     *      Does make use ofz-buffering.
     *
     * @param figures:
     *      The figures to be projected.
     * @param projectionScreenDistance:
     *      The distance of the projections screen from the eye point.
     *      It is 1.0 by default.
     */
    L2D::Lines2DZ projectFiguresZ(const L3D::Figures3D& figures, double projectionScreenDistance);


    /*
     * Description:
     *      The result is an interpolation of zaInv and zbInv,
     *      according to pxNr in relation to pxCount.
     *      The resulting value of the interpolation will then range
     *      from 1/zbInv when pxNr is 0 to 1/zaInv when
     *      pxNr is pxCount - 1.
     *
     * @param pxNr:
     *      Which of the pxCount pixels we want to find
     *      the 1/z value of.
     *      pxNr should range from 0 to pxCount - 1.
     * @param pxCount:
     *      The total number of pixels in the line.
     *      pxCount should be greater than 0.
     * @param zaInv:
     *      The value selected when pxNr is pxCount - 1.
     * @param zbInv:
     *      The value selected when pxNr is 0.
     */
    double interpolate(double pxNr, double pxCount, double zaInv, double zbInv);

    /*
     * Description:
     *      Draw a line onto the passed img::EasyImage from p0 to p1
     *      in the passed colour. z-buffering will be used to determine
     *      whether any pixel on the line should be drawn in the passed colour.
     *
     * @param x0:
     *      The x-coordinate of the first point of the to draw line.
     * @param y0:
     *      The y-coordinate of the first point of the to draw line.
     * @param z0:
     *      The z-coordinate of the first point of the to draw line.
     * @param x1:
     *      The x-coordinate of the second point of the to draw line.
     * @param y1:
     *      The y-coordinate of the second point of the to draw line.
     * @param z1:
     *      The z-coordinate of the second point of the to draw line.
     * @param color:
     *      The color which to colour the pixels of the to draw line in.
     * @param img:
     *      The img::EasImage onto which to draw the line.
     * @param ZBuffer:
     *      The z-buffer which to use to determine whether a given
     *      pixel should be recoloured.
     */
    void draw_zbuff_line(unsigned int x0, unsigned int y0, double z0,
                         unsigned int x1, unsigned int y1, double z1,
                         img::Color color, img::EasyImage& img,
                         L3D::ZBuffer& ZBuffer);

    void draw_zbuff_triag(img::EasyImage& img, L3D::ZBuffer& ZBuffer,
                          Vector3D const& A,
                          Vector3D const& B,
                          Vector3D const& C,
                          ImageSpecifications& specs,
                          L3D::LightCaster& lightCaster);

    /*
     * Description:
     *      Draw a list of projected L2D::Line2D onto a img::EasyImage and return the resulting image.
     *
     * @param lines:
     *      The list of projected lines which to draw on the generated img::EasyImage.
     * @param size:
     *      The length of the longest side of the generated img::EasyImage.
     * @param bgColor:
     *      The default color of all pixels of the generated img::EasyImage.
     *      In other words, the background colour.
     */
    img::EasyImage drawLines2D(const L2D::Lines2D& lines, double size, img::Color bgColor);

    /*
     * Description:
     *      Draw a list of projected L2D::Line2D onto a img::EasyImage and return the resulting image.
     *      z-buffering will be applied to ensure that objects closest to the eye-point
     *      will be drawn "in front of" the objects behind them.
     *
     * @param lines:
     *      The list of projected lines which to draw on the generated img::EasyImage.
     * @param size:
     *      The length of the longest side of the generated img::EasyImage.
     * @param bgColor:
     *      The default color of all pixels of the generated img::EasyImage.
     *      In other words, the background colour.
     */
    img::EasyImage drawLines2DZ(const L2D::Lines2DZ& lines, double size, img::Color bgColor);

    /*
     * Description:
     *      Draw a list of figures onto a img::EasyImage. The passed figures need to be
     *      triangulated L3D::Figure's. z-buffering will be applied to ensure that
     *      objects closest to the eye-point will be drawn "in front of" the objects behind them.
     *
     * @param figures:
     *      The list of L3D::Figure's which to draw on the generated img::EasyImage.
     * @param lightCaster:
     *      A class that will handle the lighting calculations, such that
     *      a final colour can be assigned to each drawn pixel.
     * @param size:
     *      The length of the longest side of the generated img::EasyImage.
     * @param bgColor:
     *      The default color of all pixels of the generated img::EasyImage.
     *      In other words, the background colour.
     */
    img::EasyImage drawLines2DZT(const L3D::Figures3D& figures, L3D::LightCaster& lightCaster,
                                 double size, img::Color& bgColor);


    /*
     * Description:
     *      Find the x coordinate of the intersection between the
     *      line AB and the straight  yi.
     */
    double findIntersectionX(double yi, const L2D::Point2D& A, const L2D::Point2D& B);

    /*
     * Description:
     *      Generate a imgWidth x height img::EasyImage.
     *      If neither dimensional component is 0,
     *      a red cross will be attempted tp be drawn onto the
     *      resulting image.
     *
     * @param imgWidth:
     *      The width of the generated img::EasyImage.
     * @param imgHeight:
     *      The height of the generated img::EasyImage.
     */
    img::EasyImage generateRejectionImage(unsigned int imgWidth, unsigned int imgHeight);

};


#endif //ENGINE_IMAGEGENERATOR_H
