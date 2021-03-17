//
// Created by Thomas Gueutal on 14/03/2021.
//

#ifndef ENGINE_L3D_H
#define ENGINE_L3D_H

#include "vector/vector3d.h"
#include "L2D.h"


namespace L3D {

    /*
     * Description:
     *      Create and return a 4x4 matrix that can scale a 3D point.
     *
     * @param scaleFactor:
     *      The factor by which to scale x, y and z.
     */
    Matrix scalingMatrix(double scaleFactor);

    /*
     * Description:
     *      Create and return a 4x4 matrix that rotates the
     *      y and z axes around the x axis counterclockwise.
     *
     * @param rotationAngle:
     *      The angle by which to rotate the y and z axes
     *      around the x axis counterclockwise.
     *      The input should be in radians.
     */
    Matrix rotationXMatrix(double rotationAngle);

    /*
     * Description:
     *      Create and return a 4x4 matrix that rotates the
     *      x and z axes around the y axis counterclockwise.
     *
     * @param rotationAngle:
     *      The angle by which to rotate the x and z axes
     *      around the y axis counterclockwise.
     *      The input should be in radians.
     */
    Matrix rotationYMatrix(double rotationAngle);

    /*
     * Description:
     *      Create and return a 4x4 matrix that rotates the
     *      x and y axes around the z axis counterclockwise.
     *
     * @param rotationAngle:
     *      The angle by which to rotate the x and y axes
     *      around the z axis counterclockwise.
     *      The input should be in radians.
     */
    Matrix rotationZMatrix(double rotationAngle);

    /*
     * Description:
     *      Create and return a 4x4 matrix that translates
     *      a point according to a, b and c.
     *
     * @param a:
     *      The offset added to the x-coordinate of the to translate point.
     * @param b:
     *      The offset added to the y-coordinate of the to translate point.
     * @param c:
     *      The offset added to the z-coordinate of the to translate point.
     */
    Matrix translationMatrix(double a, double b, double c);

    /*
     * Description:
     *      Create and return a 4x4 matrix that transforms
     *      a point's coordinates into the eyePoint coordinates
     *      of the passed eye's coordinate system.
     *
     * @param eye:
     *      The eys whose coordinate system we want to transform a point's
     *      coordinates into.
     */
    Matrix eyePointTransMatrix(const Vector3D& eye);

    /*
     * Description:
     *      Project a 3D point onto the plane y = -z and
     *      return that projection.
     *
     * @param point3D:
     *      The point to be projected.
     */
    inline L2D::Point2D projectPoint3D(const Vector3D& point3D, double d);



    class Figure;
    typedef std::list<Figure> Figures3D;

    /*
     * Description:
     *      This class describes a single face of a 3D figure.
     *      It consists of a list of point indexes, which refer
     *      to the points in the "points" list of the Figure instance
     *      of which this face is a part.
     *
     * @member point_indexes:
     *      The list of points that make up the face.
     *      Each point is connected to the next in the list,
     *      such that each pair forms an edge of the face.
     *      The list is circular.
     */
    class Face {

        public:

            /*
             * Description:
             *      Convert the face to a list of 2D lines,
             *      which can be drawn onto an img::EasyImage.
             *      The points of the L3D::Figure should already
             *      be in eyePoint coordinates.
             *
             * @param lines2D:
             *      The list to which to add every 2D line.
             * @param points:
             *      The list of points that the point_indexes
             *      refer to.
             * @param color:
             *      The color of each of the generated lines
             * @param projectionScreenDistance:
             *      The z distance of the projection screen from
             *      The eyePoint (the current origin).
             */
            void addToLines2D(L2D::Lines2D& lines2D,
                           const std::vector<Vector3D>& points,
                           const L2D::Color& color,
                           double projectionScreenDistance) const;

            /*
             * Description:
             *      Convert the Face to a string representation of the form
             *      ([line list], color), where a line is of the form
             *      (p1, p2).
             *
             * @param points:
             *      The list of points to which the point indexes refer.
             */
            std::string toString(const std::vector<Vector3D>& points) const;

        public:

            std::vector<int> point_indexes;
    };


    /*
     * Description:
     *      A figure is a collection of 3D points in space,
     *      where each point is part of at least one face.
     *      The points and faces are kept in lists, and the faces
     *      refer to the points list via indexes.
     *
     * @member points:
     *      The list of 3D points that make up the Figure.
     *      We make use of the provided vector3d library
     *      to represent the points as vectors.
     * @member faces:
     *      The list of faces which make up the Figure.
     *
     * @member color:
     *      The color of each of the lines of the face.
     */
    class Figure {

        public:

            explicit Figure(L2D::Color color);

            /*
             * Description:
             *      Apply a transformation to every single 3D point that is
             *      part of the L3D::Figure's points vector.
             *
             * @param M:
             *      The matrix to be applied.
             */
            void applyTransformation(const Matrix& M);

            /*
             * Description:
             *      Convert the L3D::Figure to a list of projected 2D lines.
             *      We project the 3D lines onto a projection screen/plane.
             *      The screen/plane is located at a distance of projectionScreenDistance
             *      from the eyePoint and the connecting line between the
             *      eyePoint and the origin from before the eyePoint transformation
             *      was applied is perpendicular to that plane.
             *
             * @param projectionScreenDistance:
             *      The distance to the eyePoint the projection screen is located.
             */
            L2D::Lines2D toLines2D(double projectionScreenDistance) const;

            /*
             * Description:
             *      Convert the faces that make up the face into strings
             *      and list them.
             */
            std::string toString() const;

            /*
             * Description:
             *      Find the 3D coordinates of the center of the
             *      bounding sphere of the L3D::Figure.
             */
            Vector3D boundingSphereCenter() const;

            /*
             * Description:
             *      Create and return a translation matrix that
             *      moves the L3D::Figure's bounding circle's center
             *      to the desired 3D point.
             *
             * @param targetCenter:
             *      The 3D point to which we want the bounding circle's center
             *      to move when applying the translation matrix.
             */
            Matrix centeringMatrix(const Vector3D& targetCenter) const;

            std::ostream & operator <<(std::ostream& output_stream) const;

        public:

            std::vector<Vector3D> points;
            std::vector<Face> faces;
            L2D::Color color;
    };



};


#endif //ENGINE_L3D_H
