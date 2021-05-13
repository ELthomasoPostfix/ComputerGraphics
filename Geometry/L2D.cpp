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

L2D::Point2D::Point2D() {
    this->x = 0;
    this->y = 0;
}

void L2D::Point2D::operator *=(double factor) {
    this->x *= factor;
    this->y *= factor;
}

void L2D::Point2D::operator +=(L2D::Point2D& p2) {
    this->x += p2.x;
    this->y += p2.y;
}

L2D::Point2D L2D::Point2D::operator+(const L2D::Point2D& p2) {
    return L2D::Point2D(this->x + p2.x, this->y + p2.y);
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

double L2D::Line2D::slope() const {
    if (p1.x == p2.x) {
        if (p1.y < p2.y)
            return std::numeric_limits<double>::infinity();
        return -std::numeric_limits<double>::infinity();
    }
    return (double)(p2.y - p1.y) / (double)(p2.x - p1.x);
}

void L2D::Line2D::operator *=(double factor) {
    p1 *= factor;
    p2 *= factor;
}

void L2D::Line2D::operator +=(L2D::Point2D& moveVector) {
    p1 += moveVector;
    p2 += moveVector;
}

std::string L2D::Line2D::toString() const {
    return "( (" + std::to_string(p1.x) + ", " + std::to_string(p1.y) + "), ("
               + std::to_string(p2.x) + ", " + std::to_string(p2.y) + "), ("
               + std::to_string(color.red) + ", " + std::to_string(color.green) + ", " + std::to_string(color.blue) + ") )";
}



/// ############## Line2DZ ############## ///

L2D::Line2DZ::Line2DZ(const L2D::Point2D& p1, const double z1,
                      const L2D::Point2D& p2, const double z2,
                      const double red, const double green, const double blue)
                      : L2D::Line2D(p1, p2, red, green, blue) {
    p1Z = z1;
    p2Z = z2;
}

std::string L2D::Line2DZ::toString() const {
    return "( (" + std::to_string(p1.x) + ", " + std::to_string(p1.y) + "), " + std::to_string(p1Z) + ", ("
                 + std::to_string(p2.x) + ", " + std::to_string(p2.y) + "), " + std::to_string(p2Z) + ", ("
                 + std::to_string(color.red) + ", " + std::to_string(color.green) + ", " + std::to_string(color.blue) + ") )";
}




/// ############## LGenerator ############## ///


L2D::LSystem::State::State(const double initialX, const double initialY, const double initialAngle)
                          : currentLoc(initialX, initialY)
{
    currentAngle = initialAngle;
}


/// ############## LGenerator ############## ///


L2D::LSystem::LGenerator::LGenerator() : _p1(0, 0), _p2(0, 0), _state(0, 0, 0) {}

L2D::Lines2D L2D::LSystem::LGenerator::generateLines(const ini::Configuration &configuration,
                                                       const LParser::LSystem2D& lSystem) {


    const std::vector<double> lnColor = configuration["2DLSystem"]["color"];
    img::Color lncolor(lnColor.at(0)*255.0, lnColor.at(1)*255.0, lnColor.at(2)*255.0);

    const double size = configuration["General"]["size"];

    // generate the list of lines that make up the image

    Lines2D lines;

    generateLines(lSystem, lncolor, lines);

    return lines;
}

void L2D::LSystem::LGenerator::generateLines(const LParser::LSystem2D& lSystem,
                                             img::Color& lncolor, L2D::Lines2D& lines) {

    _state.currentLoc.x = 0.0;
    _state.currentLoc.y = 0.0;
    _state.currentAngle = lSystem.get_starting_angle();
    _angle = lSystem.get_angle();

    _p1 = L2D::Point2D(0.0, 0.0);
    _p2 = L2D::Point2D(0.0, 0.0);

    // The stack should be empty after each lines generation,
    // but we empty it to avoid accidental large memory consumption
    if (!_savedStates.empty())
        _savedStates = {};

    for (const char c : lSystem.get_initiator())
    {
        recurse(c, lSystem, lSystem.get_nr_iterations(), lncolor, lines);
    }

}

void L2D::LSystem::LGenerator::recurse(const char replacedChar,
                                       const LParser::LSystem2D &lSystem,
                                       unsigned int iterations,
                                       img::Color &lncolor, L2D::Lines2D &lines) {

    // change the angle as needed
    if (replacedChar == '+') {
        _state.currentAngle += _angle;
    } else if (replacedChar == '-') {
        _state.currentAngle -= _angle;
    } else if (replacedChar == '(') {
        _savedStates.push(_state);
    } else if (replacedChar == ')') {
        if (_savedStates.empty())
            std::cerr << "During the generation of a 2D LSystem, a ')' was encountered."
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
        _p1 = L2D::Point2D(_state.currentLoc.x, _state.currentLoc.y);

        // adjust the current coordinates
        _state.currentLoc.x += std::cos(toRadians(_state.currentAngle));
        _state.currentLoc.y += std::sin(toRadians(_state.currentAngle));

        // end point of the new line
        _p2 = L2D::Point2D(_state.currentLoc.x, _state.currentLoc.y);

        // add a new line to the list
        lines.emplace_back(L2D::Line2D(_p1, _p2, lncolor.red, lncolor.green, lncolor.blue));

    } else {

        // adjust the current coordinates
         _state.currentLoc.x += std::cos(_state.currentAngle);
         _state.currentLoc.y += std::sin(_state.currentAngle);

    }


}
