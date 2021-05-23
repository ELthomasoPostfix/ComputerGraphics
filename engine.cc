#include "utils/ini_configuration.h"
#include "utils/l_parser.h"
#include "ImageGeneration/ImageGenerator.h"



img::EasyImage generate_image(const ini::Configuration &configuration, const std::string& iniFilePath) {
    const std::string& type = configuration["General"]["type"].as_string_or_die();

    if (type == "2DLSystem") {

        // retrieve configuration attributes
        double size = configuration["General"]["size"].as_double_or_die();
        std::vector<double> bgColor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        img::Color bgcolor(bgColor.at(0)*255.0, bgColor.at(1)*255.0, bgColor.at(2)*255.0);

        // create the l_parser and read the LSystem2D input file into it

        LParser::LSystem2D lSystem2D;
        const std::string& LFile = configuration["2DLSystem"]["inputfile"].as_string_or_die();
        const std::string L2DFile = truncateFileName(iniFilePath) + LFile;
        // std::cout << LFile << std::endl; // TODO delete ???
        // std::cout << "\t" << L2DFile << std::endl; // TODO delete ???
        std::ifstream inputStream(L2DFile);
        inputStream >> lSystem2D;
        inputStream.close();

        L2D::LSystem::LGenerator lSystemGenerator;
        L2D::Lines2D lines = lSystemGenerator.generateLines(configuration, lSystem2D);
        return ImageGenerator::drawLines2D(lines, size, bgcolor);
    }
    else if (type == "Wireframe" || type == "ZBufferedWireframe" || type == "ZBuffering" || type == "LightedZBuffering") {

        // Whether or not triangulation should be performed on the L3D::Figures.
        const bool triangulate = type == "ZBuffering" || type == "LightedZBuffering";

        // retrieve configuration attributes
        double size = configuration["General"]["size"].as_double_or_die();

        std::vector<double> bgColor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        img::Color bgcolor(bgColor.at(0)*255.0, bgColor.at(1)*255.0, bgColor.at(2)*255.0);

        L3D::LightCaster lightCaster = ImageGenerator::parseLights(configuration);

        // Parse and create the list of figures.
        // These figures are composed of Vector3D points.
        L3D::Figures3D figures = ImageGenerator::parseFigures(configuration, iniFilePath, triangulate);

        if (figures.empty())
            goto reject;

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

        lightCaster.applyTransformation(eyeTrans);
        lightCaster.eye = eye;

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



        if (type == "ZBuffering" || type == "LightedZBuffering") {
            return ImageGenerator::drawLines2DZT(figures, lightCaster, size, bgcolor);

        } else if (type == "ZBufferedWireframe") {
            // generate the 2D lines
            L2D::Lines2DZ lines2D = ImageGenerator::projectFiguresZ(figures, 1.0);

            return ImageGenerator::drawLines2DZ(lines2D, size, bgcolor);
        } else if (type == "Wireframe") {
            // generate the 2D lines
            L2D::Lines2D lines2D = ImageGenerator::projectFigures(figures, 1.0);

            return ImageGenerator::drawLines2D(lines2D, size, bgcolor);
        }
    }

    reject:
    return ImageGenerator::generateRejectionImage(1024, 1024);
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
