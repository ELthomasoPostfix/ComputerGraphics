//
// Created by Thomas Gueutal on 10/05/2021.
//

#include "ImageGenerator.h"


L3D::Figures3D ImageGenerator::parseFigures(const ini::Configuration& configuration,
                                            const std::string& iniFilePath,
                                            const bool triangulate) {

    L3D::Figures3D figures = {};

    for (unsigned int figureIndex = 0; figureIndex < (unsigned int) configuration["General"]["nrFigures"].as_int_or_die(); figureIndex++) {

        // setup configuration variables
        const std::string figureName = "Figure" + std::to_string(figureIndex);
        const std::string type = configuration[figureName]["type"].as_string_or_die();

        std::vector<double> Color = configuration[figureName]["color"].as_double_tuple_or_die();
        L2D::Color color(Color.at(0)*255.0, Color.at(1)*255.0, Color.at(2)*255.0);

        if (type == "LineDrawing") {
            figures.emplace_back(L3D::Figure::createLineDrawingFigure(color, configuration, figureName));

        } else if (type == "3DLSystem") {

            LParser::LSystem3D lSystem3D;
            const std::string& LFile = configuration[figureName]["inputfile"].as_string_or_die();
            const std::string L3DFile = truncateFileName(iniFilePath) + LFile;
            std::cout << "\t" << L3DFile << std::endl;
            std::ifstream inputStream(L3DFile);
            inputStream >> lSystem3D;
            inputStream.close();

            L3D::LSystem::LGenerator lSystemGenerator;
            figures.emplace_back(lSystemGenerator.generateFigure(configuration, color, lSystem3D));

        } else if (type == "Cube") {
            figures.emplace_back(L3D::Figure::createCube(color, triangulate));
        } else if (type == "Tetrahedron") {
            figures.emplace_back(L3D::Figure::createTetrahedron(color));
        } else if (type == "Octahedron") {
            figures.emplace_back(L3D::Figure::createOctahedron(color));
        } else if (type == "Icosahedron") {
            figures.emplace_back(L3D::Figure::createIcosahedron(color));
        } else if (type == "Dodecahedron") {
            figures.emplace_back(L3D::Figure::createDodecahedron(color, triangulate));
        } else if (type == "Cylinder") {
            figures.emplace_back(L3D::Figure::createCylinder(color, configuration, figureName, triangulate));
        } else if (type == "Cone") {
            figures.emplace_back(L3D::Figure::createCone(color, configuration, figureName, triangulate));
        } else if (type == "Sphere") {
            figures.emplace_back(L3D::Figure::createSphere(color, configuration, figureName));
        } else if (type == "Torus") {
            figures.emplace_back(L3D::Figure::createTorus(color, configuration, figureName, triangulate));
        }

    }

    return figures;
}



L2D::Lines2D ImageGenerator::projectFigures(const L3D::Figures3D& figures, const double projectionScreenDistance) {

    L2D::Lines2D lines2D = {};

    for (const L3D::Figure& fig : figures) {
        L2D::Lines2D newLines = fig.toLines2D(projectionScreenDistance);
        lines2D.insert(lines2D.end(), newLines.begin(), newLines.end());
    }

    return lines2D;
}

L2D::Lines2DZ ImageGenerator::projectFiguresZ(const L3D::Figures3D& figures, const double projectionScreenDistance) {

    L2D::Lines2DZ lines2DZ = {};

    for (const L3D::Figure& fig : figures) {
        L2D::Lines2DZ newLines = fig.toLines2DZ(projectionScreenDistance);
        lines2DZ.insert(lines2DZ.end(), newLines.begin(), newLines.end());
    }

    return lines2DZ;
}



double ImageGenerator::interpolate(const double pxNr, const double pxCount, const double zaInv, const double zbInv) {
    // TODO is this correct ?????
    return ((pxNr / (pxCount - 1)) * zaInv) + ((1 - (pxNr / (pxCount - 1))) * zbInv);
}

void ImageGenerator::draw_zbuff_line(unsigned int x0, unsigned int y0, const double z0,
                                     unsigned int x1, unsigned int y1, const double z1,
                                     img::Color color, img::EasyImage& img,
                                     L3D::ZBuffer& ZBuffer) {

    assert(x0 < img.get_width() && y0 < img.get_height());
    assert(x1 < img.get_width() && y1 < img.get_height());

    double z0Inv = 1.0 / z0;
    double z1Inv = 1.0 / z1;
    unsigned int pxCount;     // nr of pixels to draw
    unsigned int pxNr;
    double zInv;

    if (x0 == x1)
    {
        pxCount = std::max(y0, y1) - std::min(y0, y1) + 1;
        pxNr = pxCount - 1;

        // The pixels will be coloured in the order of the smallest to the largest of y0, y1,
        // so we make sure that the smallest of the two's zInv is in z0Inv.
        if (y1 < y0)
            std::swap(z0Inv, z1Inv);

        // special case for x0 == x1
        for (unsigned int i = std::min(y0, y1); i <= std::max(y0, y1); i++)
        {
            zInv = interpolate(pxNr, pxCount, z0Inv, z1Inv);
            if (ZBuffer.replace(x0, i, zInv))
                img(x0, i) = color;
            pxNr--;
        }
    }
    else if (y0 == y1)
    {
        pxCount = std::max(x0, x1) - std::min(x0, x1) + 1;
        pxNr = pxCount - 1;

        // The pixels will be coloured in the order of the smallest to the largest of x0, x1,
        // so we make sure that the smallest of the two's zInv is in z0Inv.
        if (x1 < x0)
            std::swap(z0Inv, z1Inv);

        // special case for y0 == y1
        for (unsigned int i = std::min(x0, x1); i <= std::max(x0, x1); i++)
        {
            zInv = interpolate(pxNr, pxCount, z0Inv, z1Inv);
            if (ZBuffer.replace(i, y0, zInv))
                img(i, y0) = color;
            pxNr--;
        }
    }
    else
    {
        unsigned int x;
        unsigned int y;

        if (x0 > x1)
        {
            // flip points if x1>x0: we want x0 to have the lowest value
            std::swap(x0, x1);
            std::swap(y0, y1);
            std::swap(z0Inv, z1Inv);
        }
        double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
        if (-1.0 <= m && m <= 1.0)
        {
            pxCount = x1 - x0 + 1;
            pxNr = pxCount - 1;

            for (unsigned int i = 0; i <= (x1 - x0); i++)
            {
                x = x0 + i;
                y = (unsigned int) roundToInt((double) y0 + m * (double) i);

                zInv = interpolate(pxNr, pxCount, z0Inv, z1Inv);
                if (ZBuffer.replace(x, y, zInv))
                    img(x, y) = color;
                pxNr--;
            }
        }
        else if (m > 1.0)
        {
            pxCount = y1 - y0 + 1;
            pxNr = pxCount - 1;

            for (unsigned int i = 0; i <= (y1 - y0); i++)
            {
                x = (unsigned int) roundToInt((double) x0 + ((double) i / m));
                y = y0 + i;

                zInv = interpolate(pxNr, pxCount, z0Inv, z1Inv);
                if (ZBuffer.replace(x, y, zInv))
                    img(x, y) = color;
                pxNr--;
            }
        }
        else if (m < -1.0)
        {
            pxCount = y0 - y1 + 1;
            pxNr = pxCount - 1;

            for (unsigned int i = 0; i <= (y0 - y1); i++)
            {
                x = (unsigned int) roundToInt((double) x0 - ((double) i / m));
                y = y0 - i;

                zInv = interpolate(pxNr, pxCount, z0Inv, z1Inv);
                if (ZBuffer.replace(x, y, zInv))
                    img(x, y) = color;
                pxNr--;
            }
        }
    }
}

void ImageGenerator::draw_zbuff_triag(img::Color color, img::EasyImage& img, L3D::ZBuffer& ZBuffer,
                                      Vector3D const& A,
                                      Vector3D const& B,
                                      Vector3D const& C,
                                      ImageSpecifications& specs) {

    Vector3D u = B - A;
    Vector3D v = C - A;
    Vector3D w = Vector3D();
    w.x = u.y*v.z - u.z*v.y;
    w.y = u.z*v.x - u.x*v.z;
    w.z = u.x*v.y - u.y*v.x;

    const double k = w.x*A.x + w.y*A.y + w.z*A.z;
    const double dzdx = w.x / ( - specs.projectionScreenDistance * k);
    const double dzdy = w.y / ( - specs.projectionScreenDistance * k);

    L2D::Point2D moveVector = specs.getMoveVector();
    const double infty = std::numeric_limits<double>::infinity();

    L2D::Point2D AProj = L3D::projectPoint3D(A, specs.projectionScreenDistance) + moveVector;
    L2D::Point2D BProj = L3D::projectPoint3D(B, specs.projectionScreenDistance) + moveVector;
    L2D::Point2D CProj = L3D::projectPoint3D(C, specs.projectionScreenDistance) + moveVector;

    const unsigned int yMin = roundToInt(std::min(std::min(AProj.y, BProj.y), CProj.y) + 0.5);
    const unsigned int yMax = roundToInt(std::max(std::max(AProj.y, BProj.y), CProj.y) - 0.5);

    const double xg = (AProj.x + BProj.x + CProj.x) / 3.0;
    const double yg = (AProj.y + BProj.y + CProj.y) / 3.0;
    const double zgInv = 1.0/3.0 * (1.0/A.z + 1.0/B.z + 1.0/C.z);

    std::vector<L2D::Point2D*> edgePoints = {
            &AProj, &BProj,
            &AProj, &CProj,
            &BProj, &CProj
    };

    std::vector<double> intersections = {
              infty,   infty,   infty,
            - infty, - infty, - infty
    };

    const L2D::Point2D* p1;
    const L2D::Point2D* p2;

    for (unsigned int yCurr = yMin; yCurr <= yMax; yCurr++) {
        // Find the two intersections.
        for (unsigned int i = 0; i < 5; i += 2) {
            p1 = edgePoints[i];
            p2 = edgePoints[i + 1];

            if (p1 != nullptr) {
                if (p1->y != p2->y) {
                    if ((yCurr - p1->y)*(yCurr - p2->y) <= 0) {
                        const double intersec = findIntersectionX(yCurr, *p1, *p2);
                        intersections[i]   = intersec;
                        intersections[i+3] = intersec;
                    }
                } else      // If [p1 p2] is a horizontal edge, then we ignore it when finding any intersections.
                    edgePoints[i] = nullptr;
            }
        }

        const unsigned int xL = (unsigned int) std::min(std::min(intersections[0], intersections[1]), intersections[2]);
        const unsigned int xR = (unsigned int) std::max(std::max(intersections[3], intersections[4]), intersections[5]);
        intersections[0] = infty;
        intersections[1] = infty;
        intersections[2] = infty;
        intersections[3] = - infty;
        intersections[4] = - infty;
        intersections[5] = - infty;

        const double zInvCte = 1.0001 * zgInv;

        for (unsigned int xCurr = xL; xCurr <= xR; xCurr++) {
            const double zInv = zInvCte + (((double) xCurr) - xg) * dzdx + (((double) yCurr) - yg) * dzdy;
            // try coloring z-buffer[xCurr][yCurr] with zInv = zInv and color color
            if (ZBuffer.replace(xCurr, yCurr, zInv)) {
                // img.draw_line(xCurr, yCurr, xCurr, yCurr, color);
                img(xCurr, yCurr) = color;
            }
        }
    }
}


img::EasyImage ImageGenerator::drawLines2D(const L2D::Lines2D& lines, const double size, img::Color bgColor) {

    ImageSpecifications specs = ImageSpecifications(lines, size);

    L2D::Point2D moveVector = specs.getMoveVector();

    img::EasyImage image(roundToInt(specs.imageX), roundToInt(specs.imageY), bgColor);

    // draw the finalized lines
    for (const L2D::Line2D& line : lines)
    {
        L2D::Line2D lineCopy = line;

        // scale the line by a factor of d
        lineCopy *= specs.projectionScreenDistance;

        // move the line to the correct location
        lineCopy += moveVector;

        // draw the finalized line
        image.draw_line(roundToInt(lineCopy.p1.x), roundToInt(lineCopy.p1.y),
                        roundToInt(lineCopy.p2.x), roundToInt(lineCopy.p2.y),
                        img::Color(roundToInt(lineCopy.color.red), roundToInt(lineCopy.color.green), roundToInt(lineCopy.color.blue)));
    }

    return image;
}

img::EasyImage ImageGenerator::drawLines2DZ(const L2D::Lines2DZ& lines, const double size, img::Color bgColor) {

    ImageSpecifications specs = ImageSpecifications(lines, size);

    L2D::Point2D moveVector = specs.getMoveVector();

    img::EasyImage image(roundToInt(specs.imageX), roundToInt(specs.imageY), bgColor);

    L3D::ZBuffer ZBuffer = L3D::ZBuffer(image.get_width(), image.get_height());

    // draw the finalized lines
    for (const L2D::Line2DZ& line : lines)
    {
        L2D::Line2DZ lineCopy = line;

        // scale the line by a factor of d
        lineCopy *= specs.projectionScreenDistance;

        // move the line to the correct location
        lineCopy += moveVector;

        // draw the finalized line
        draw_zbuff_line(roundToInt(lineCopy.p1.x), roundToInt(lineCopy.p1.y), lineCopy.p1Z,
                        roundToInt(lineCopy.p2.x), roundToInt(lineCopy.p2.y), lineCopy.p2Z,
                        img::Color(roundToInt(lineCopy.color.red), roundToInt(lineCopy.color.green), roundToInt(lineCopy.color.blue)),
                        image,
                        ZBuffer);
    }

    return image;
}

img::EasyImage ImageGenerator::drawLines2DZT(const L3D::Figures3D& figures, const double size, img::Color& bgColor) {

    ImageSpecifications specs = ImageSpecifications(figures, size);

    img::EasyImage image(roundToInt(specs.imageX), roundToInt(specs.imageY), bgColor);

    L3D::ZBuffer ZBuffer = L3D::ZBuffer(image.get_width(), image.get_height());

    for (const L3D::Figure& fig : figures) {
        for (const L3D::Face& face : fig.faces) {
            draw_zbuff_triag(img::Color(roundToInt(fig.color.red), roundToInt(fig.color.green), roundToInt(fig.color.blue)),
                             image,
                             ZBuffer,
                             fig.points.at(face.point_indexes.at(0)),
                             fig.points.at(face.point_indexes.at(1)),
                             fig.points.at(face.point_indexes.at(2)),
                             specs);
        }
    }

    return image;
}



double ImageGenerator::findIntersectionX(const double yi, const L2D::Point2D& A, const L2D::Point2D& B) {

    // TODO  if A = B = C then draw a line

    double xLeft = std::min(A.x, B.x);
    double xRight = std::max(A.x, B.x);
    double yUpper = std::min(A.y, B.y);
    double yLower = std::max(A.y, B.y);

    double slope;
    if (A.x == B.x)
        slope = - std::numeric_limits<double>::infinity();      // Distinguishing between +infinity and -infinity
    else                                                        // is not helpful. Both cases result in intersection
        slope = (B.y - A.y) / (B.x - A.x);                      // at {xLeft, yi}.

    if (slope == std::numeric_limits<double>::infinity())
        return xLeft;
    else if (slope < 0)
        return xLeft + ((yUpper - yi) * (xRight - xLeft) / (yUpper - yLower));
    else
        return xLeft + ((yi - yLower) * (xRight - xLeft) / (yUpper - yLower));
}
