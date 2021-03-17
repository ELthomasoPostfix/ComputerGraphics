//
// Created by Thomas Gueutal on 22/02/2021.
//

#ifndef ENGINE_L2D_H
#define ENGINE_L2D_H

#include "utils/easy_image.h"
#include "utils/ini_configuration.h"
#include "utils/l_parser.h"
#include "utils/utils.h"

#include <cmath>
#include <list>
#include <stack>

namespace L2D {

    class Line2D;
    using Lines2D = std::list<Line2D>;


    class Color {
        public:
            Color(double red, double green, double blue);

        public:
            double red;
            double green;
            double blue;
    };


    class Point2D {
        public:
            Point2D(double x, double y);

            void operator *=(double factor);

            void operator +=(L2D::Point2D& p2);

    public:
            double x;
            double y;
    };


    class Line2D {
        public:
            Line2D(const Point2D& p1, const Point2D& p2,
                   double red, double green, double blue);

            Line2D(double x1, double y1, double x2, double y2,
                   double red, double green, double blue);

            void operator *=(double factor);

            void operator +=(L2D::Point2D& moveVector);

            std::string toString() const;

        public:
            Point2D p1;
            Point2D p2;
            L2D::Color color;
    };


    namespace LSystem {


        /*
         * Description:
         *       This struct represents the two mutable, base components of the LGenerator.
         *       They support the generation process.
         *
         * @member currentLoc:
         *      The current x-coordinate and y-coordinate in 2D space.
         * @member currentAngle:
         *      The current angle under which to to move/draw.
         *      !! Stored in degrees !!
         */
        struct State {
            public:

                State(double initialX, double initialY, double initialAngle);

            public:
                L2D::Point2D currentLoc;
                double currentAngle;
        };



        /*
         *
         * Description:
         *      The LGenerator class takes in a ini::Configuration and a LParser::LSystem2D.
         *      The configuration is used to construct the returned img::EasyImage,
         *      and to draw the generated lines onto that image. The lsystem is used
         *      to generate the list of lines to be drawn onto the image.
         *
         * @member _state:
         *      The current location and the current angle, represented as the state of the LGenerator.
         *      The state contains the two mutable, base attributes that support the generation process.
         * @member angle:
         *      The angle by which to modify the current angle.
         *      !! Stored in degrees !!
         * @member _p1:
         *      The start point of the to generate line.
         *      _p1 is (re)assigned just before adjusting the current coordinates.
         * @member _p2:
         *      The end point of the to generate line.
         *      _p2 is (re)assigned just after adjusting the current coordinates.
         * @member _savedStates:
         *      It is the set of states saved during the generation of the lines.
         *      They are only valid within the generation they were saved in.
         *      If a '[' is encountered, then _state will be saved on the stack.
         *      When the corresponding ']' is encountered, then the current state
         *      will be set to the last saved state.
         */

        class LGenerator {
        private:

            /*
             * Description:
             *      Use the initiator of the provided lSystem to generate the lines
             *      to be drawn.
             *
             * @param lSystem:
             *      The lSystem based on which to generate the lines.
             * @param lncolor:
             *      The color every line is drawn in.
             * @param lines:
             *      The list to which to add the lines.
             */
            void generateLines(const LParser::LSystem2D& lSystem,
                               img::Color& lncolor, Lines2D& lines);


            /*
             * Description:
             *      Based on the replacedChar (not '+', '-', '(', ')'), if the nr of iterations left is > 0,
             *      then recurse onto the replacement rule of replacedChar.
             *      Else, adjust the angle ('+', '-'), save or load the state ('(', ')') or
             *      possibly generate a line.
             *
             * @param replacedChar:
             *      The character on which we are recursing.
             * @param lSystem:
             *      The lSystem we retrieve the replacement rules from.
             * @param iterations
             *      The nr of recurse iterations left.
             * @param lncolor:
             *      The color which to draw the line in.
             * @param lines:
             *      The list to which to add the generated lines.
             */
            void recurse(char replacedChar,
                         const LParser::LSystem2D& lSystem,
                         unsigned int iterations,
                         img::Color& lncolor, Lines2D& lines);

        public:

            LGenerator();


            /*
             * Description:
             *      Generate an img::EasyImage based on the passed
             *      ini::Configuration image requirements and lSystem.
             *
             * @param configuration:
             *      The requirements for the to generate img::EasyImage
             * @param lSystem:
             *      The lSystem based on which to generate the lines to draw
             *      onto the img::EasyImage.
             */
            img::EasyImage generateImage(const ini::Configuration &configuration, const LParser::LSystem2D& lSystem);


        private:

            State _state;

            double _angle;

            L2D::Point2D _p1;
            L2D::Point2D _p2;

            std::stack<State> _savedStates;
        };


    }

}
#endif //ENGINE_L2D_H
