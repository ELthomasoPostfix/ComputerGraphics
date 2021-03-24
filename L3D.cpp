//
// Created by Thomas Gueutal on 14/03/2021.
//

#include "L3D.h"


Matrix L3D::scalingMatrix(const double scaleFactor) {

    Matrix scaleMatrix;

    // assign the scale factor to the first three diagonal elements
    scaleMatrix(1, 1) = scaleFactor;
    scaleMatrix(2, 2) = scaleFactor;
    scaleMatrix(3, 3) = scaleFactor;

    return scaleMatrix;
}

Matrix L3D::rotationXMatrix(const double rotationAngle) {

    Matrix rotationMatrix;

    // construct the rotation matrix
    rotationMatrix(2, 2) =   std::cos(rotationAngle);
    rotationMatrix(2, 3) =   std::sin(rotationAngle);

    rotationMatrix(3, 2) = - std::sin(rotationAngle);
    rotationMatrix(3, 3) =   std::cos(rotationAngle);

    return rotationMatrix;
}

Matrix L3D::rotationYMatrix(const double rotationAngle) {

    Matrix rotationMatrix;

    // construct the rotation matrix
    rotationMatrix(1, 1) =   std::cos(rotationAngle);
    rotationMatrix(1, 3) = - std::sin(rotationAngle);

    rotationMatrix(3, 1) =   std::sin(rotationAngle);
    rotationMatrix(3, 3) =   std::cos(rotationAngle);

    return rotationMatrix;
}

Matrix L3D::rotationZMatrix(const double rotationAngle) {

    Matrix rotationMatrix;

    // construct the rotation matrix
    rotationMatrix(1, 1) =   std::cos(rotationAngle);
    rotationMatrix(1, 2) =   std::sin(rotationAngle);

    rotationMatrix(2, 1) = - std::sin(rotationAngle);
    rotationMatrix(2, 2) =   std::cos(rotationAngle);

    return rotationMatrix;
}

Matrix L3D::translationMatrix(const double a, const double b, const double c) {

    Matrix translationMatrix;

    // assign the offsets to the first three columns of the last row
    translationMatrix(4, 1) = a;
    translationMatrix(4, 2) = b;
    translationMatrix(4, 3) = c;

    return translationMatrix;
}

Matrix L3D::eyePointTransMatrix(const Vector3D& eye) {

    // convert the cartesian eye coordinates to polar coordinates
    const double r     = std::sqrt(std::pow(eye.x, 2) + std::pow(eye.y, 2) + std::pow(eye.z, 2));
    const double theta = std::atan2(eye.y, eye.x);
    const double phi   = std::acos(eye.z / r);

    // create the eye transformation matrix
    Matrix eyeTrans = L3D::rotationZMatrix(theta + (M_PI / 2.0))
                    * L3D::rotationXMatrix(phi)
                    * L3D::translationMatrix(0, 0, -r);

    Matrix eyeTransManual;
    eyeTransManual(1, 1) = - std::sin(theta);
    eyeTransManual(1, 2) = - std::cos(theta) * std::cos(phi);
    eyeTransManual(1, 3) =   std::cos(theta) * std::sin(phi);

    eyeTransManual(2, 1) =   std::cos(theta);
    eyeTransManual(2, 2) = - std::sin(theta) * std::cos(phi);
    eyeTransManual(2, 3) =   std::sin(theta) * std::sin(phi);

    eyeTransManual(3, 2) =   std::sin(phi);
    eyeTransManual(3, 3) =   std::cos(phi);

    eyeTransManual(4, 3) =   -r;

    return eyeTransManual;
}

inline L2D::Point2D L3D::projectPoint3D(const Vector3D& point3D, const double d) {

    return { d * point3D.x / (-point3D.z),
             d * point3D.y / (-point3D.z) };
}





void L3D::Face::addToLines2D(L2D::Lines2D& lines2D,
                             const std::vector<Vector3D> &points,
                             const L2D::Color& color,
                             const double projectionScreenDistance) const {

    // if the face consists of only two points (a singe line),
    // then the list should not be seen as circular
    if (point_indexes.size() == 2) {
        const Vector3D &p1 = points.at(point_indexes.at(0));
        const Vector3D &p2 = points.at(point_indexes.at(1));

        lines2D.emplace_back(L2D::Line2D(projectPoint3D(p1, projectionScreenDistance),
                                         projectPoint3D(p2, projectionScreenDistance),
                                         color.red, color.green, color.blue));
    } else {
        // TODO if a point index pair is found in multiple faces, then we make a line for each occurrence !!! =! good
        // TODO ==> add a condition: toDrawDuplicate
        unsigned i2 = 1;
        for (unsigned int i1 = 0; i1 < point_indexes.size(); i1++) {

            // retrieve the 3D points that need to be converted
            const Vector3D &p1 = points.at(point_indexes.at(i1));
            const Vector3D &p2 = points.at(point_indexes.at(i2));

            // project a 3D line of two 3D points onto a 2D line
            lines2D.emplace_back(L2D::Line2D(projectPoint3D(p1, projectionScreenDistance),
                                             projectPoint3D(p2, projectionScreenDistance),
                                             color.red, color.green, color.blue));

            // i2 needs to be able to wrap around, as the indexes list is circular
            i2 = (i2 == point_indexes.size() - 1) ? 0 : i2 + 1;
        }
    }
}

std::string L3D::Face::toString(const std::vector<Vector3D>& points) const {
    std::string res;

    // add the lines
    if (point_indexes.size() == 2) {
        Vector3D p1 = points.at(point_indexes.at(0));
        Vector3D p2 = points.at(point_indexes.at(1));

        res += "\t( (" + std::to_string(p1.x) + ", " + std::to_string(p1.y) + ", " + std::to_string(p1.z) + "), ("
                     + std::to_string(p2.x) + ", " + std::to_string(p2.y) + ", " + std::to_string(p2.z) + ") )\n";

    } else {
        unsigned int i2 = 1;
        for (unsigned int i1 = 0; i1 < point_indexes.size(); i1++) {

            // retrieve the 3D points that need to be converted
            const Vector3D &p1 = points.at(point_indexes.at(i1));
            const Vector3D &p2 = points.at(point_indexes.at(i2));

            // add the line to the result string
            res += "\t( (" + std::to_string(p1.x) + ", " + std::to_string(p1.y) + ", " + std::to_string(p1.z) + "), ("
                         + std::to_string(p2.x) + ", " + std::to_string(p2.y) + ", " + std::to_string(p2.z) + ") )\n";

            // i2 needs to be able to wrap around, as the indexes list is circular
            i2 = (i2 == point_indexes.size() - 1) ? 0 : i2 + 1;
        }
    }

    return res;
}



L3D::Figure::Figure(L2D::Color color) : color(color) {}

void L3D::Figure::applyTransformation(const Matrix &M) {

    for (Vector3D& point : points) {
        point *= M;
    }
}

L2D::Lines2D L3D::Figure::toLines2D(const double projectionScreenDistance) const {

    L2D::Lines2D lines2D = {};

    for (const L3D::Face& face : faces) {
        face.addToLines2D(lines2D, points, color, projectionScreenDistance);
    }

    return lines2D;
}

std::string L3D::Figure::toString() const {

    std::string res;

    for (unsigned int index = 0; index < faces.size(); ++index) {
        res += "Face" + std::to_string(index) + "(\n";
        res += "\t" + faces.at(index).toString(points);
        res += ")\n";
    }

    return res;
}

std::ostream &L3D::Figure::operator<<(std::ostream &output_stream) const {

    output_stream << this->toString();

    return output_stream;
}



