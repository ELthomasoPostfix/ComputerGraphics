#include "utils/ini_configuration.h"
#include "utils/l_parser.h"
#include "L3D.h"


inline double interpolate(const double pxNr, const double pxCount, const double zaInv, const double zbInv) {
    return ((pxNr / (pxCount - 1)) * zaInv) + ((1 - (pxNr / (pxCount - 1))) * zbInv);
}

void draw_zbuff_line(unsigned int x0, unsigned int y0, const double z0,
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
            std::swap(z0Inv, z1Inv);     // TODO: commenting makes z_buffered_wireframes118.ini correct
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

L3D::Figures3D parseFigures(const ini::Configuration& configuration, const std::string& iniFilePath = "") {

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
            figures.emplace_back(L3D::Figure::createCube(color));
        } else if (type == "Tetrahedron") {
            figures.emplace_back(L3D::Figure::createTetrahedron(color));
        } else if (type == "Octahedron") {
            figures.emplace_back(L3D::Figure::createOctahedron(color));
        } else if (type == "Icosahedron") {
            figures.emplace_back(L3D::Figure::createIcosahedron(color));
        } else if (type == "Dodecahedron") {
            figures.emplace_back(L3D::Figure::createDodecahedron(color));
        } else if (type == "Cylinder") {
            figures.emplace_back(L3D::Figure::createCylinder(color, configuration, figureName));
        } else if (type == "Cone") {
            figures.emplace_back(L3D::Figure::createCone(color, configuration, figureName));
        } else if (type == "Sphere") {
            figures.emplace_back(L3D::Figure::createSphere(color, configuration, figureName));
        } else if (type == "Torus") {
            figures.emplace_back(L3D::Figure::createTorus(color, configuration, figureName));
        }

    }

    return figures;
}

L2D::Lines2D projectFigures(const L3D::Figures3D& figures, const double projectionScreenDistance) {

    L2D::Lines2D lines2D = {};

    for (const L3D::Figure& fig : figures) {
        L2D::Lines2D newLines = fig.toLines2D(projectionScreenDistance);
        lines2D.insert(lines2D.end(), newLines.begin(), newLines.end());
    }

    return lines2D;
}

L2D::Lines2DZ projectFiguresZ(const L3D::Figures3D& figures, const double projectionScreenDistance) {

    L2D::Lines2DZ lines2DZ = {};

    for (const L3D::Figure& fig : figures) {
        L2D::Lines2DZ newLines = fig.toLines2DZ(projectionScreenDistance);
        lines2DZ.insert(lines2DZ.end(), newLines.begin(), newLines.end());
    }

    return lines2DZ;
}

img::EasyImage drawLines2D(const L2D::Lines2D& lines, const double size, img::Color bgColor) {

    double minX =  std::numeric_limits<double>::infinity();
    double minY =  std::numeric_limits<double>::infinity();
    double maxX = -std::numeric_limits<double>::infinity();
    double maxY = -std::numeric_limits<double>::infinity();

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
    double imageX = size * (rangeX / std::max(rangeX, rangeY));

    double imageY = size * (rangeY / std::max(rangeX, rangeY));

    // find the factor by which to scale every line
    double d = 0.95 * (imageX / rangeX);

    // Calculate the rescaled center of the set of lines
    double DCx = d * ((minX + maxX) / 2.0);
    double DCy = d * ((minY + maxY) / 2.0);

    // Find the new center of the lines
    double dx = (imageX / 2.0) - DCx;
    double dy = (imageY / 2.0) - DCy;

    L2D::Point2D moveVector(dx, dy);

    img::EasyImage image(roundToInt(imageX), roundToInt(imageY), bgColor);



    // draw the finalized lines
    for (const L2D::Line2D& line : lines)
    {
        L2D::Line2D lineCopy = line;

        // scale the line by a factor of d
        lineCopy *= d;

        // move the line to the correct location
        lineCopy += moveVector;

        // draw the finalized line
        image.draw_line(roundToInt(lineCopy.p1.x), roundToInt(lineCopy.p1.y),
                        roundToInt(lineCopy.p2.x), roundToInt(lineCopy.p2.y),
                  img::Color(roundToInt(lineCopy.color.red), roundToInt(lineCopy.color.green), roundToInt(lineCopy.color.blue)));
    }

    return image;
}

img::EasyImage drawLines2DZ(const L2D::Lines2DZ& lines, const double size, img::Color bgColor) {

    double minX =  std::numeric_limits<double>::infinity();
    double minY =  std::numeric_limits<double>::infinity();
    double maxX = -std::numeric_limits<double>::infinity();
    double maxY = -std::numeric_limits<double>::infinity();

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
    double imageX = size * (rangeX / std::max(rangeX, rangeY));

    double imageY = size * (rangeY / std::max(rangeX, rangeY));

    // find the factor by which to scale every line
    double d = 0.95 * (imageX / rangeX);

    // Calculate the rescaled center of the set of lines
    double DCx = d * ((minX + maxX) / 2.0);
    double DCy = d * ((minY + maxY) / 2.0);

    // Find the new center of the lines
    double dx = (imageX / 2.0) - DCx;
    double dy = (imageY / 2.0) - DCy;

    L2D::Point2D moveVector(dx, dy);

    img::EasyImage image(roundToInt(imageX), roundToInt(imageY), bgColor);

    L3D::ZBuffer ZBuffer = L3D::ZBuffer(image.get_width(), image.get_height());

    // draw the finalized lines
    for (const L2D::Line2DZ& line : lines)
    {
        L2D::Line2DZ lineCopy = line;

        // scale the line by a factor of d
        lineCopy *= d;

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

img::EasyImage generate_image(const ini::Configuration &configuration, const std::string& iniFilePath)
{
    const std::string& type = configuration["General"]["type"].as_string_or_die();

    if (type == "2DLSystem") {

        // create the l_parser and read the LSystem2D input file into it

        LParser::LSystem2D lSystem2D;
        const std::string& LFile = configuration["2DLSystem"]["inputfile"].as_string_or_die();
        const std::string L2DFile = truncateFileName(iniFilePath) + LFile;
        std::cout << LFile << std::endl; // TODO delete ???
        std::cout << "\t" << L2DFile << std::endl; // TODO delete ???
        std::ifstream inputStream(L2DFile);
        inputStream >> lSystem2D;
        inputStream.close();

        L2D::LSystem::LGenerator lSystemGenerator;
        return lSystemGenerator.generateImage(configuration, lSystem2D);
    }
    else if (type == "Wireframe" || type == "ZBufferedWireframe") {

        // retrieve configuration attributes
        double size = configuration["General"]["size"].as_double_or_die();

        std::vector<double> bgColor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        img::Color bgcolor(bgColor.at(0)*255.0, bgColor.at(0)*255.0, bgColor.at(0)*255.0);

        // Parse and create the list of figures.
        // These figures are composed of Vector3D points.
        L3D::Figures3D figures = parseFigures(configuration, iniFilePath);

        // apply any needed transformations
        Matrix fullTrans;       // The complete transformation that should be applied to all points

        double rotateX = 0;     // Rotation angle around the x-axis (degrees).
        double rotateY = 0;     // Rotation angle around the y-axis (degrees).
        double rotateZ = 0;     // Rotation angle around the z-axis (degrees).
        std::vector<double> center = {};    // The final center location of the figure, after all transformations
        double scale = 0;       // The factor by which to scale the points
        std::string figureName; // the name by which to access the figure configuration attributes

        std::vector<double> eyeCoordinates = configuration["General"]["eye"].as_double_tuple_or_die();
        Vector3D eye = Vector3D::point(eyeCoordinates.at(0), eyeCoordinates.at(1), eyeCoordinates.at(2));
        Matrix eyeTrans = L3D::eyePointTransMatrix(eye);

        unsigned int figIndex = 0;
        for (L3D::Figure& figure : figures) {

            // retrieve configuration data
            figureName = "Figure" + std::to_string(figIndex);
            scale   = configuration[figureName]["scale"].as_double_or_die();
            rotateX = configuration[figureName]["rotateX"].as_double_or_die();
            rotateY = configuration[figureName]["rotateY"].as_double_or_die();
            rotateZ = configuration[figureName]["rotateZ"].as_double_or_die();
            center  = configuration[figureName]["center"].as_double_tuple_or_die();

            // Create all sub matrices
            Matrix S   = L3D::scalingMatrix(scale);             // scaling matrix
            Matrix RX  = L3D::rotationXMatrix(toRadians(rotateX));  // rotation around the x axis
            Matrix RY  = L3D::rotationYMatrix(toRadians(rotateY));  // rotation around the y axis
            Matrix RZ  = L3D::rotationZMatrix(toRadians(rotateZ));  // rotation around the z axis
            Matrix C   = L3D::translationMatrix(center.at(0), center.at(1), center.at(2));  // the intended center for the figure

            // assemble and apply the full transformation matrix
            fullTrans = S * RX * RY * RZ * C * eyeTrans;

            figure.applyTransformation(fullTrans);

            figIndex++;
        }

        if (type == "ZBufferedWireframe") {
            // generate the 2D lines
            L2D::Lines2DZ lines2D = projectFiguresZ(figures, 1.0);

            return drawLines2DZ(lines2D, size, bgcolor);
        } else {
            // generate the 2D lines
            L2D::Lines2D lines2D = projectFigures(figures, 1.0);

            return drawLines2D(lines2D, size, bgcolor);
        }
    }


    return img::EasyImage();
}

int main(int argc, char const* argv[])
{
    int retVal;
        try
        {
                for(int i = 1; i < argc; ++i)
                {
                        ini::Configuration conf;
                        try
                        {
                            std::cerr << argv[i] << std::endl;
                                std::ifstream fin(argv[i]);
                                fin >> conf;
                                fin.close();

                        }
                        catch(ini::ParseException& ex)
                        {
                                std::cerr << "Error parsing file: " << argv[i] << ": " << ex.what() << std::endl;
                                retVal = 1;
                                continue;
                        }

                        img::EasyImage image = generate_image(conf, argv[i]);
                        if(image.get_height() > 0 && image.get_width() > 0)
                        {
                                std::string fileName(argv[i]);
                                std::string::size_type pos = fileName.rfind('.');
                                if(pos == std::string::npos)
                                {
                                        //filename does not contain a '.' --> append a '.bmp' suffix
                                        fileName += ".bmp";
                                }
                                else
                                {
                                        fileName = fileName.substr(0,pos) + ".bmp";
                                }
                                try
                                {
                                        std::ofstream f_out(fileName.c_str(),std::ios::trunc | std::ios::out | std::ios::binary);
                                        f_out << image;

                                }
                                catch(std::exception& ex)
                                {
                                        std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                                        retVal = 1;
                                }
                        }
                        else
                        {
                                std::cout << "Could not generate image for " << argv[i] << std::endl;
                        }
                }
        }
        catch(const std::bad_alloc &exception)
        {
    		//When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
    		//Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
    		//(Unless of course you are already consuming the maximum allowed amount of memory)
    		//If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
		//mark the test as failed while in reality it just needed a bit more memory
                std::cerr << "Error: insufficient memory" << std::endl;
                retVal = 100;
        }
        return retVal;
}
