//
// Created by Thomas Gueutal on 22/02/2021.
//

#include "L2D.h"



/// ############## Color ############## ///

L2D::Color::Color(const double red, const double green, const double blue) {
    this->red = red;
    this->green = green;
    this->blue = blue;
}



/// ############## Point2D ############## ///

L2D::Point2D::Point2D(const double x, const double y) {
    this->x = x;
    this->y = y;
}

void L2D::Point2D::operator *=(double factor) {
    this->x *= factor;
    this->y *= factor;
}

void L2D::Point2D::operator +=(L2D::Point2D& p2) {
    this->x += p2.x;
    this->y += p2.y;
}



/// ############## Line2D ############## ///

L2D::Line2D::Line2D(const Point2D& p1, const Point2D& p2,
                    double red, double green, double blue) :
    p1(p1.x, p1.y),
    p2(p2.x, p2.y),
    color(red, green, blue) {}

L2D::Line2D::Line2D(double x1, double y1, double x2, double y2,
                    double red, double green, double blue) :
    p1(x1, y1),
    p2(x2, y2),
    color(red, green, blue) {}

void L2D::Line2D::operator *=(double factor) {
    p1 *= factor;
    p2 *= factor;
}

void L2D::Line2D::operator +=(L2D::Point2D& moveVector) {
    p1 += moveVector;
    p2 += moveVector;
}



/// ############## LGenerator ############## ///


L2D::LSystem::State::State(const double initialX, const double initialY, const double initialAngle)
                          : _currentLoc(initialX, initialY)
{
    _currentAngle = initialAngle;
}


/// ############## LGenerator ############## ///


L2D::LSystem::LGenerator::LGenerator() : _p1(0, 0), _p2(0, 0), _state(0, 0, 0) {}

// function found in engine.cc
// creates an img::EasyImage and draws the passed list of lines onto it
img::EasyImage drawLines2D(const L2D::Lines2D& lines, double size, img::Color bgColor);


img::EasyImage L2D::LSystem::LGenerator::generateImage(const ini::Configuration &configuration,
                                                       const LParser::LSystem2D& lSystem) {

    // properties for creating the easy image

    const std::vector<double> bgColor = configuration["General"]["backgroundcolor"];
    img::Color bgcolor(bgColor.at(0)*255.0, bgColor.at(1)*255.0, bgColor.at(2)*255.0);

    const std::vector<double> lnColor = configuration["2DLSystem"]["color"];
    img::Color lncolor(lnColor.at(0)*255.0, lnColor.at(1)*255.0, lnColor.at(2)*255.0);

    const double size = configuration["General"]["size"];

    // generate the list of lines that make up the image

    Lines2D lines;

    generateLines(lSystem, lncolor, lines);

    return drawLines2D(lines, size, bgcolor);
}

void L2D::LSystem::LGenerator::generateLines(const LParser::LSystem2D& lSystem,
                                             img::Color& lncolor, Lines2D& lines) {

    _state._currentLoc.x = 0.0;
    _state._currentLoc.y = 0.0;
    _state._currentAngle = lSystem.get_starting_angle();
    _angle = lSystem.get_angle();

    _p1 = L2D::Point2D(0.0, 0.0);
    _p2 = L2D::Point2D(0.0, 0.0);

    // The stack should be empty after each lines generation,
    // but we empty it to avoid accidental large memory consumption
    if (!_savedStates.empty())
        _savedStates = {};

    for (const char c : lSystem.get_initiator())
    {
        std::cout << c << std::endl;
        recurse(c, lSystem, lSystem.get_nr_iterations(), lncolor, lines);
    }

}

void L2D::LSystem::LGenerator::recurse(const char replacedChar,
                                       const LParser::LSystem2D &lSystem,
                                       unsigned int iterations,
                                       img::Color &lncolor, L2D::Lines2D &lines) {

    // change the angle as needed
    if (replacedChar == '+') {
        _state._currentAngle += _angle;
    } else if (replacedChar == '-') {
        _state._currentAngle -= _angle;
    } else if (replacedChar == '(') {
        _savedStates.push(_state);
    } else if (replacedChar == ')') {
        if (_savedStates.empty())
            std::cerr << "During the generation of a LSystem, a ']' was encountered."
                      << "There was, however, no saved stated to reload to." << std::endl;

        _state = _savedStates.top();
        _savedStates.pop();
    }
    // if recursion necessary, do it
    else if (iterations > 0)
    {
        for (const char c : lSystem.get_replacement(replacedChar))
        {
            recurse(c, lSystem, iterations-1, lncolor, lines);
        }
    }
    // no more recursion needed, so start adding lines
    else if (lSystem.draw(replacedChar)) {

        // starting point of the new line
        _p1 = L2D::Point2D(_state._currentLoc.x, _state._currentLoc.y);

        // adjust the current coordinates
        _state._currentLoc.x += std::cos(toRadians(_state._currentAngle));
        _state._currentLoc.y += std::sin(toRadians(_state._currentAngle));

        // end point of the new line
        _p2 = L2D::Point2D(_state._currentLoc.x, _state._currentLoc.y);

        // add a new line to the list
        lines.emplace_back(L2D::Line2D(_p1, _p2, lncolor.red, lncolor.green, lncolor.blue));

    } else {

        // adjust the current coordinates
         _state._currentLoc.x += std::cos(toRadians(_state._currentAngle));
         _state._currentLoc.y += std::sin(toRadians(_state._currentAngle));

    }


}

