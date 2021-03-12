#include "utils/ini_configuration.h"
#include "utils/l_parser.h"
#include "L2D.h"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <assert.h>



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

    std::cout << "TYPE : " << type << std::endl;

    if (type == "2DLSystem") {

        // create the l_parser and read the LSystem2D input file into it

        LParser::LSystem2D lSystem2D;
        const std::string& LFile = configuration["2DLSystem"]["inputfile"].as_string_or_die();
        std::cout << LFile << std::endl; // TODO
        std::ifstream inputStream("../build/l_systems/" + LFile);
        inputStream >> lSystem2D;
        inputStream.close();

        L2D::LSystem::LGenerator lSystemGenerator;
        return lSystemGenerator.generateImage(configuration, lSystem2D);

        // return drawLines2D(lines, size, bgcolor);
    }
    else {


        double size = 500;



        // TODO  test

        img::Color bgcolor;
        bgcolor.red = 2;
        bgcolor.green = 0;
        bgcolor.blue = 0;
        L2D::Color lnColor(0, 255, 0);

        // TODO  test
        std::cout << "bgColor = (" << (unsigned int)bgcolor.red << ", " << (unsigned int)bgcolor.green << ", " << (unsigned int)bgcolor.blue << ")" << std::endl;
        std::cout << "lnColor = (" << lnColor.red << ", " << lnColor.green << ", " << lnColor.blue << ")" << std::endl;
        // TODO  test

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
