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
        bool addedFigure = false;

        L2D::Color ambient = L2D::Color::black();
        L2D::Color diffuse = L2D::Color::black();
        L2D::Color specular = L2D::Color::black();
        double reflectionCoefficient = 0;

        // Retrieve the refraction components if they should be present
        // in the ini::Configuration.
        if (configuration["General"]["type"].as_string_or_die() == "LightedZBuffering") {
            std::vector<double> ar = configuration[figureName]["ambientReflection"].as_double_tuple_or_default({0.0, 0.0, 0.0});
            ambient = L2D::Color::colorClamp(ar[0], ar[1], ar[2]);
            std::vector<double> dr = configuration[figureName]["diffuseReflection"].as_double_tuple_or_default({0.0, 0.0, 0.0});
            diffuse = L2D::Color::colorClamp(dr[0], dr[1], dr[2]);
            std::vector<double> sr = configuration[figureName]["specularReflection"].as_double_tuple_or_default({0.0, 0.0, 0.0});
            specular = L2D::Color::colorClamp(sr[0], sr[1], sr[2]);

            configuration[figureName]["reflectionCoefficient"].as_double_if_exists(reflectionCoefficient);
        }
        // Else treat the otherwise present color as the ambient component.
        else {
            std::vector<double> color = configuration[figureName]["color"].as_double_tuple_or_die();
            ambient = L2D::Color::colorClamp(color[0], color[1], color[2]);

        }

        if (type == "LineDrawing") {
            figures.emplace_back(L3D::Figure::createLineDrawingFigure(configuration, figureName));
            addedFigure = true;
        } else if (type == "3DLSystem") {

            LParser::LSystem3D lSystem3D;
            const std::string& LFile = configuration[figureName]["inputfile"].as_string_or_die();
            const std::string L3DFile = truncateFileName(iniFilePath) + LFile;
            // std::cout << "\t" << L3DFile << std::endl;   // TODO  print
            std::ifstream inputStream(L3DFile);
            inputStream >> lSystem3D;
            inputStream.close();

            L3D::LSystem::LGenerator lSystemGenerator;
            figures.emplace_back(lSystemGenerator.generateFigure(configuration, lSystem3D));
            addedFigure = true;
        } else if (type == "Cube") {
            figures.emplace_back(L3D::Figure::createCube(triangulate));
            addedFigure = true;
        } else if (type == "Tetrahedron") {
            figures.emplace_back(L3D::Figure::createTetrahedron());
            addedFigure = true;
        } else if (type == "Octahedron") {
            figures.emplace_back(L3D::Figure::createOctahedron());
            addedFigure = true;
        } else if (type == "Icosahedron") {
            figures.emplace_back(L3D::Figure::createIcosahedron());
            addedFigure = true;
        } else if (type == "Dodecahedron") {
            figures.emplace_back(L3D::Figure::createDodecahedron(triangulate));
            addedFigure = true;
        } else if (type == "Cylinder") {
            figures.emplace_back(L3D::Figure::createCylinder(configuration, figureName, triangulate));
            addedFigure = true;
        } else if (type == "Cone") {
            figures.emplace_back(L3D::Figure::createCone(configuration, figureName, triangulate));
            addedFigure = true;
        } else if (type == "Sphere") {
            figures.emplace_back(L3D::Figure::createSphere(configuration, figureName));
            addedFigure = true;
        } else if (type == "Torus") {
            figures.emplace_back(L3D::Figure::createTorus(configuration, figureName, triangulate));
            addedFigure = true;
        }

        if (addedFigure) {
            L3D::Figure& justAdded = figures.back();
            justAdded.ambientReflectivity  = ambient;
            justAdded.diffuseReflectivity  = diffuse;
            justAdded.specularReflectivity = specular;
            justAdded.reflectionCoefficient = reflectionCoefficient;
        }

        addedFigure = false;

    }

    return figures;
}


L3D::LightCaster ImageGenerator::parseLights(const ini::Configuration& configuration) {

    L3D::LightCaster lightCaster;
    std::string type = configuration["General"]["type"];

    if (type == "LightedZBuffering") {

        unsigned int lightCount = configuration["General"]["nrLights"].as_int_or_die();

        for (unsigned int lightIndex = 0; lightIndex < lightCount; lightIndex++) {

            const std::string lightName = "Light" + std::to_string(lightIndex);

            std::vector<double> ai = configuration[lightName]["ambientLight"].as_double_tuple_or_default({0.0, 0.0, 0.0});
            L2D::Color ambientIntensity = L2D::Color::colorClamp(ai.at(0), ai.at(1), ai.at(2));
            std::vector<double> di = configuration[lightName]["diffuseLight"].as_double_tuple_or_default({0.0, 0.0, 0.0});
            L2D::Color diffuseIntensity = L2D::Color::colorClamp(di.at(0), di.at(1), di.at(2));
            std::vector<double> si = configuration[lightName]["specularLight"].as_double_tuple_or_default({0.0, 0.0, 0.0});
            L2D::Color specularIntensity = L2D::Color::colorClamp(si.at(0), si.at(1), si.at(2));

            bool infinity = configuration[lightName]["infinity"].as_bool_or_default(false);

            std::vector<double> vector;
            const std::string vectorName = infinity ? "direction" : "location";

            if (configuration[lightName][vectorName].as_double_tuple_if_exists(vector)) {
                if (infinity) {
                    L3D::InfLight* light = new L3D::InfLight(
                            ambientIntensity, diffuseIntensity, specularIntensity,
                            Vector3D::vector(vector[0], vector[1], vector[2]));
                    lightCaster.lightSources.emplace_back(light);
                } else {

                    double spotAngle;
                    if (configuration[lightName]["spotAngle"].as_double_if_exists(spotAngle)) {
                        L3D::PointLight* pointLight = new L3D::PointLight(
                                ambientIntensity, diffuseIntensity, specularIntensity,
                                Vector3D::point(vector[0], vector[1], vector[2]),
                                toRadians(clamp(spotAngle, 0.0, 90.0)));
                        lightCaster.lightSources.emplace_back(pointLight);
                    } else {
                        L3D::PointLight* pointLight = new L3D::PointLight(
                                ambientIntensity, diffuseIntensity, specularIntensity,
                                Vector3D::point(vector[0], vector[1], vector[2]),
                                -1.0);
                        lightCaster.lightSources.emplace_back(pointLight);
                    }
                }
            } else {
                L3D::Light* ambientLight = new L3D::Light(ambientIntensity, diffuseIntensity, specularIntensity);
                lightCaster.lightSources.emplace_back(ambientLight);
            }

        }

    } else {
        // Generate the default ambient light for non lighting configurations.
        L2D::Color ambientIntensity = L2D::Color::colorClamp(1.0, 1.0, 1.0);
        L2D::Color diffuseIntensity = L2D::Color::black();
        L2D::Color specularIntensity = L2D::Color::black();
        L3D::Light* ambientLight = new L3D::Light(ambientIntensity, diffuseIntensity, specularIntensity);

        lightCaster.lightSources.emplace_back(ambientLight);
    }

    return lightCaster;
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

        // The pixels will be coloured in the order of the smallest to the largest of {y0, y1},
        // so we make sure that the smallest of the two y's zInv is in z0Inv.
        if (y1 < y0)
            std::swap(z0Inv, z1Inv);

        // special case for x0 == x1
        for (unsigned int i = std::min(y0, y1); i <= std::max(y0, y1); i++)
        {
            // interpolate chooses z0Inv if pxNr = pxCount - 1
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
            // interpolate chooses z0Inv if pxNr = pxCount - 1
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

                // interpolate chooses z0Inv if pxNr = pxCount - 1
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

                // interpolate chooses z0Inv if pxNr = pxCount - 1
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

                // interpolate chooses z0Inv if pxNr = pxCount - 1
                zInv = interpolate(pxNr, pxCount, z0Inv, z1Inv);
                if (ZBuffer.replace(x, y, zInv))
                    img(x, y) = color;
                pxNr--;
            }
        }
    }
}

void ImageGenerator::draw_zbuff_triag(img::EasyImage& img, L3D::ZBuffer& ZBuffer,
                                      Vector3D const& A,
                                      Vector3D const& B,
                                      Vector3D const& C,
                                      ImageSpecifications& specs,
                                      L3D::LightCaster& lightCaster) {

    Vector3D triangleSideB = B - A;
    Vector3D triangleSideC = C - A;
    Vector3D normal = Vector3D::cross(triangleSideB, triangleSideC);

    Vector3D normalized = Vector3D::normalise(normal);
    lightCaster.recalculateInfResults(normalized);

    const double k = Vector3D::dot(normal, A);
    const double dzdx = normal.x / ( - specs.projectionScreenDistance * k);
    const double dzdy = normal.y / ( - specs.projectionScreenDistance * k);

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


    L2D::Point2D* edgePoints[6] = {
            &AProj, &BProj,
            &AProj, &CProj,
            &BProj, &CProj
    };

    const L2D::Point2D* p1;
    const L2D::Point2D* p2;


    for (unsigned int yCurr = yMin; yCurr <= yMax; yCurr++) {

        double intersections[6] = {
                infty,   infty,   infty,
                - infty, - infty, - infty
        };

        // Find the two intersections.
        for (unsigned int i = 0; i < 5; i += 2) {
            p1 = edgePoints[i];
            p2 = edgePoints[i + 1];

            if (p1 != nullptr) {
                if (p1->y != p2->y) {
                    if ((yCurr - p1->y)*(yCurr - p2->y) <= 0) {
                        const double intersec = findIntersectionX(yCurr, *p1, *p2);
                        unsigned int index1 = i/2;
                        unsigned int index2 = index1 + 3;
                        intersections[index1] = intersec; // 0 1 2
                        intersections[index2] = intersec; // 3 4 5
                    }
                } else      // If [p1 p2] is a horizontal edge, then we ignore it when finding any intersections.
                    edgePoints[i] = nullptr;
            }
        }

        const unsigned int xL = roundToInt(std::min(std::min(intersections[0], intersections[1]), intersections[2]) + 0.5);
        const unsigned int xR = roundToInt(std::max(std::max(intersections[3], intersections[4]), intersections[5]) - 0.5);

        const double zInvCte = 1.0001 * zgInv;

        for (unsigned int xCurr = xL; xCurr <= xR; xCurr++) {
            const double zInv = zInvCte + (((double) xCurr) - xg) * dzdx + (((double) yCurr) - yg) * dzdy;
            // try coloring z-buffer[xCurr][yCurr] with zInv = zInv and color color
            if (ZBuffer.replace(xCurr, yCurr, zInv)) {

                lightCaster.recalculatePointResults(normalized,
                                                (double) xCurr - specs.dx,
                                                (double) yCurr - specs.dy,
                                                1.0/zInv, specs.projectionScreenDistance);
                img(xCurr, yCurr) = lightCaster.getClampedResultColor().toImageColor();
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
                        lineCopy.color.toImageColor());
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
                        lineCopy.color.toImageColor(),
                        image,
                        ZBuffer);
    }

    return image;
}

img::EasyImage ImageGenerator::drawLines2DZT(const L3D::Figures3D& figures, L3D::LightCaster& lightCaster,
                                             const double size, img::Color& bgColor) {

    ImageSpecifications specs = ImageSpecifications(figures, size);

    img::EasyImage image(roundToInt(specs.imageX), roundToInt(specs.imageY), bgColor);

    L3D::ZBuffer ZBuffer = L3D::ZBuffer(image.get_width(), image.get_height());
    for (const L3D::Figure& fig : figures) {
        lightCaster.recalculateAmbientResult(fig);
        lightCaster.setReflectivityComponents(fig);
        for (const L3D::Face& face : fig.faces) {
            draw_zbuff_triag(image,
                             ZBuffer,
                             fig.points.at(face.point_indexes.at(0)),
                             fig.points.at(face.point_indexes.at(1)),
                             fig.points.at(face.point_indexes.at(2)),
                             specs,
                             lightCaster);
        }
    }

    return image;
}



double ImageGenerator::findIntersectionX(const double yi, const L2D::Point2D& A, const L2D::Point2D& B) {

    // TODO  if A = B = C then draw a line

    double xLeft = std::min(A.x, B.x);
    double xRight = std::max(A.x, B.x);
    double yUpper = std::max(A.y, B.y);
    double yLower = std::min(A.y, B.y);

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


img::EasyImage ImageGenerator::generateRejectionImage(const unsigned int imgWidth, const unsigned int imgHeight) {
    img::EasyImage rejImg = img::EasyImage(imgWidth, imgHeight, L2D::Color::black().toImageColor());
    unsigned int xRight = imgWidth-1, yUpper = imgHeight-1;

    if (imgWidth == 0 || imgHeight == 0)
        return rejImg;

    rejImg.draw_line(0, 0, xRight, yUpper, L2D::Color(1.0, 0.0, 0.0).toImageColor());
    rejImg.draw_line(0, yUpper, xRight, 0, L2D::Color(1.0, 0.0, 0.0).toImageColor());
    return rejImg;
}

