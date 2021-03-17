#include "utils/ini_configuration.h"
#include "utils/l_parser.h"
#include "L3D.h"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <assert.h>


L3D::Figures3D createFigures(const ini::Configuration& configuration) {

    L3D::Figures3D figures = {};

    for (unsigned int figureIndex = 0; figureIndex < (unsigned int) configuration["General"]["nrFigures"].as_int_or_die(); figureIndex++) {

        // setup configuration variables
        const std::string figureName = "Figure" + std::to_string(figureIndex);

        std::vector<double> Color = configuration[figureName]["color"].as_double_tuple_or_die();
        L2D::Color color(Color.at(0)*255.0, Color.at(1)*255.0, Color.at(2)*255.0);


        // create a new figure
        L3D::Figure newFigure = L3D::Figure(color);

        // add all the 3D points to the L3D::Figure
        std::vector<double> point = {};
        for (unsigned int pointIndex = 0; pointIndex < (unsigned int) configuration[figureName]["nrPoints"].as_int_or_die(); pointIndex++) {

            point = configuration[figureName]["point" + std::to_string(pointIndex)];

            newFigure.points.emplace_back(Vector3D::point(point.at(0), point.at(1), point.at(2)));
        }

        // add all faces to the figure
        // TODO find the actual faces, instead of making each line a face of its own
        std::vector<int> line = {};
        for (unsigned int lineIndex = 0; lineIndex < (unsigned int) configuration[figureName]["nrLines"].as_int_or_die(); lineIndex++) {

            L3D::Face newFace;

            line = configuration[figureName]["line" + std::to_string(lineIndex)].as_int_tuple_or_die();
            newFace.point_indexes.emplace_back(line.at(0));
            newFace.point_indexes.emplace_back(line.at(1));

            newFigure.faces.emplace_back(newFace);
        }

        figures.emplace_back(newFigure);
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

img::EasyImage drawLines2D(const L2D::Lines2D& lines, const double size, img::Color bgColor) {

    double minX = 32767.0;
    double minY = 32767.0;
    double maxX = -32768.0;
    double maxY = -32768.0;

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
    assert(rangeX >= 0.0 ||rangeY >= 0.0);

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

img::EasyImage generate_image(const ini::Configuration &configuration)
{
    const std::string& type = configuration["General"]["type"].as_string_or_die();

    if (type == "2DLSystem") {

        // create the l_parser and read the LSystem2D input file into it

        LParser::LSystem2D lSystem2D;
        const std::string& LFile = configuration["2DLSystem"]["inputfile"].as_string_or_die();
        std::cout << LFile << std::endl; // TODO
        std::ifstream inputStream("../ini/l_systems/" + LFile);
        inputStream >> lSystem2D;
        inputStream.close();

        L2D::LSystem::LGenerator lSystemGenerator;
        return lSystemGenerator.generateImage(configuration, lSystem2D);
    }
    else if (type == "Wireframe") {

        // retrieve configuration attributes
        double size = configuration["General"]["size"].as_double_or_die();

        std::vector<double> bgColor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        img::Color bgcolor(bgColor.at(0)*255.0, bgColor.at(0)*255.0, bgColor.at(0)*255.0);

        // create the list of figures
        L3D::Figures3D figures = createFigures(configuration);

        // apply any needed transformations
        Matrix fullTrans;       // The complete transformation that should be applied to all points

        double rotateX = 0;     // How many degrees a figure should be rotated around the x-axis.
        double rotateY = 0;     // How many degrees a figure should be rotated around the y-axis.
        double rotateZ = 0;     // How many degrees a figure should be rotated around the z-axis.
        std::vector<double> center = {};    // The final center location of the figure, after all transformations
        double scale = 0;       // The factor by which to scale the points
        std::string figureName; // the name by which to access the figure configuration attributes

        std::vector<double> eyeCoordinates = configuration["General"]["eye"].as_double_tuple_or_die();
        Vector3D eye = Vector3D::point(eyeCoordinates.at(0), eyeCoordinates.at(1), eyeCoordinates.at(2));
        Matrix eyeTrans = L3D::eyePointTransMatrix(eye);

        unsigned int figIndex = 0;
        for (L3D::Figure& figure : figures) {

            Vector3D targetCenter = Vector3D::point(0,0,0);   // center the figure on around (0, 0, 0)

            // retrieve configuration data
            figureName = "Figure" + std::to_string(figIndex);
            scale   = configuration[figureName]["scale"].as_double_or_die();
            rotateX = configuration[figureName]["rotateX"].as_double_or_die();
            rotateY = configuration[figureName]["rotateY"].as_double_or_die();
            rotateZ = configuration[figureName]["rotateZ"].as_double_or_die();
            center  = configuration[figureName]["center"].as_double_tuple_or_die();

            // Create all sub matrices
            Matrix OCM = figure.centeringMatrix(targetCenter);  // origin centering matrix
            Matrix S   = L3D::scalingMatrix(scale);             // scaling matrix
            Matrix RX  = L3D::rotationXMatrix(toRadians(rotateX));  // rotation around the x axis
            Matrix RY  = L3D::rotationYMatrix(toRadians(rotateY));  // rotation around the y axis
            Matrix RZ  = L3D::rotationZMatrix(toRadians(rotateZ));  // rotation around the z axis
            Matrix C   = L3D::translationMatrix(center.at(0), center.at(1), center.at(2));  // the intended center for the figure

            // assemble and apply the full transformation matrix
            fullTrans = OCM * S * RX * RY * RZ * C * eyeTrans;  // ALL FIGURES ROTATED SLIGHTLY
            // fullTrans = S * RX * RY * RZ * C * eyeTrans;        // NEARLY WORKS, only small differences

            figure.applyTransformation(fullTrans);

            figIndex++;
        }

        // generate the 2D lines
        L2D::Lines2D lines2D = projectFigures(figures, 1);

        return drawLines2D(lines2D, size, bgcolor);
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

                        img::EasyImage image = generate_image(conf);
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
