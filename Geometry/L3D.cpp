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

L2D::Point2D L3D::projectPoint3D(const Vector3D& point3D, const double d) {

    return { d * point3D.x / (-point3D.z),
             d * point3D.y / (-point3D.z) };
}

Vector3D L3D::recreatePoint3D(const double x, const double y, const double z, const double d) {
    return Vector3D::point(x * (-z) / d, y * (-z) / d, z);
}





L3D::Light::Light(const L2D::Color &ambient, const L2D::Color &diffuse, const L2D::Color &specular) {
    this->_ambientIntensity  = ambient;
    this->_diffuseIntensity  = diffuse;
    this->_specularIntensity = specular;
}

L2D::Color
L3D::Light::getDiffuseContribution(const Vector3D &normalizedNormal, const L2D::Color &diffuseReflectivity) const {
    return L2D::Color::black();
}

L2D::Color L3D::Light::getDiffuseContribution(const Vector3D &normalizedNormal, const L2D::Color &diffuseReflectivity,
                                              double x, double y, double z, double d) const {
    return L2D::Color::black();
}

L2D::Color
L3D::Light::getSpecularContribution(const Vector3D &normalizedNormal, const L2D::Color &diffuseReflectivity, double x,
                                    double y, double z, double d) const {
    return L2D::Color::black();
}

const L2D::Color &L3D::Light::getAmbientIntensity() const {
    return this->_ambientIntensity;
}

const L2D::Color &L3D::Light::getDiffuseIntensity() const {
    return this->_diffuseIntensity;
}

const L2D::Color &L3D::Light::getSpecularIntensity() const {
    return this->_specularIntensity;
}

unsigned int L3D::Light::getType() const {
    return 0;
}

bool L3D::Light::isInfLight() const {
    return this->getType() == 1;
}

bool L3D::Light::isPointLight() const {
    return this->getType() == 2;
}

void L3D::Light::applyTransformation(const Matrix &matrix) {}





L3D::InfLight::InfLight(const L2D::Color &ambient, const L2D::Color &diffuse, const L2D::Color &specular,
                        const Vector3D &ldVector) :
             L3D::Light(ambient, diffuse, specular) {

    this->_ldVector = ldVector;
}

L2D::Color L3D::InfLight::getDiffuseContribution(const Vector3D &normalizedNormal, const L2D::Color& diffuseReflectivity) const {
    if (!diffuseReflectivity.nonZero())
        return L2D::Color::black();

    const double cosAlpha = Vector3D::dot(normalizedNormal, - Vector3D::normalise(this->_ldVector));

    if (cosAlpha <= 0.0)
        return L2D::Color::black();

    return diffuseReflectivity * this->getDiffuseIntensity() * cosAlpha;
}

L2D::Color L3D::InfLight::getSpecularContribution(const Vector3D &normalizedNormal, const L2D::Color &specularReflectivity,
                                                  double x, double y, double z, double d, const Vector3D& eye) const {

    if (!specularReflectivity.nonZero())
        return L2D::Color::black();

    Vector3D l = - Vector3D::normalise(this->_ldVector);
    const double cosAlpha = Vector3D::dot(normalizedNormal, l);

    if (cosAlpha <= 0.0)
        return L2D::Color::black();

    Vector3D r = (normalizedNormal * 2.0 * cosAlpha) - l;
    Vector3D originalPoint = L3D::recreatePoint3D(x, y, z, d);
    Vector3D eyeConnector = eye - originalPoint;

    return specularReflectivity * this->getSpecularIntensity() * cosAlpha;
}

void L3D::InfLight::applyTransformation(const Matrix &matrix) {
    this->_ldVector *= matrix;
}

unsigned int L3D::InfLight::getType() const {
    return 1;
}






L3D::PointLight::PointLight(const L2D::Color &ambient, const L2D::Color &diffuse, const L2D::Color &specular,
                            const Vector3D &location, const double spotAngle) :
                            L3D::Light(ambient, diffuse, specular) {

    this->_location = location;
    this->_spotAngle = spotAngle;
}

L2D::Color
L3D::PointLight::getDiffuseContribution(const Vector3D &normalizedNormal, const L2D::Color &diffuseReflectivity) const {
    return L2D::Color::black();
}

L2D::Color
L3D::PointLight::getDiffuseContribution(const Vector3D &normalizedNormal, const L2D::Color &diffuseReflectivity,
                                        const double x, const double y, const double z, const double d) const {

    if (!diffuseReflectivity.nonZero())
        return L2D::Color::black();

    const bool spotAngleProvided = this->_spotAngle >= 0.0;

    Vector3D originalPoint = L3D::recreatePoint3D(x, y, z, d);
    Vector3D ldVector = this->_location - originalPoint;

    const double cosAlpha = Vector3D::dot(normalizedNormal,  Vector3D::normalise(ldVector));
    double spotAngle;

    if (spotAngleProvided) {
        spotAngle = std::cos(this->_spotAngle);
        if (cosAlpha < spotAngle)
            return L2D::Color::black();
    }
    if (cosAlpha <= 0.0)
        return L2D::Color::black();

    double angleModifier = spotAngleProvided ?
            1.0 - (1.0 - cosAlpha)/(1.0 - spotAngle) :
            cosAlpha;

    return diffuseReflectivity * this->getDiffuseIntensity() * angleModifier;

}

void L3D::PointLight::applyTransformation(const Matrix &matrix) {
    this->_location *= matrix;
}

unsigned int L3D::PointLight::getType() const {
    return 2;
}









L3D::LightCaster::LightCaster() :
        _ambientResult(0.0, 0.0, 0.0),
        _infDiffuseResult(0.0, 0.0, 0.0),
        _pointDiffuseResult(0.0, 0.0, 0.0),
        _infSpecularResult(0.0, 0.0, 0.0),
        _pointSpecularResult(0.0, 0.0, 0.0),
        _diffuseReflectivity(nullptr),
        _specularReflectivity(nullptr),
        _reflectionCoefficient(0) {}

L3D::LightCaster::~LightCaster() {
    for (Light* light : lightSources)
        delete light;
}

void L3D::LightCaster::setReflectivityComponents(const L3D::Figure &figure) {
    this->_diffuseReflectivity = &figure.diffuseReflectivity;
    this->_specularReflectivity = &figure.specularReflectivity;
    this->_reflectionCoefficient = figure.reflectionCoefficient;
}

void L3D::LightCaster::recalculateAmbientResult(const L3D::Figure& figure) {
    this->_ambientResult = L2D::Color::black();
    for (const L3D::Light* lightSource : lightSources)
        this->_ambientResult += lightSource->getAmbientIntensity();

    this->_ambientResult *= figure.ambientReflectivity;
}

void L3D::LightCaster::recalculateInfResults(const Vector3D &normalizedNormal) {
    // Reset the infinite diffuse component and add in all the contributions
    // of the L3D::InfLight light sources.
    this->_infDiffuseResult = L2D::Color::black();
    for (const L3D::Light* light : lightSources) {
        if (light->isInfLight()) {
            this->_infDiffuseResult += light->getDiffuseContribution(normalizedNormal, *this->_diffuseReflectivity);
            // this->_infSpecularResult += light->getS
        }
    }
}

void L3D::LightCaster::recalculatePointResults(const Vector3D &normalizedNormal, const double x,
                                              const double y, const double z, const double d) {
    // Reset the point diffuse component and add in all the contributions
    // of the L3D::PointLight light sources.
    this->_pointDiffuseResult = L2D::Color::black();
    for (const L3D::Light* light : lightSources) {
        if (light->isPointLight())
            this->_pointDiffuseResult += light->getDiffuseContribution(normalizedNormal, *this->_diffuseReflectivity, x, y, z, d);
    }
}

L2D::Color L3D::LightCaster::getResultColor() const {
    return this->_ambientResult + this->_infDiffuseResult + this->_pointDiffuseResult +
    this->_infSpecularResult + this->_pointSpecularResult;
}

L2D::Color L3D::LightCaster::getClampedResultColor() const {
    L2D::Color unClamped = this->_ambientResult + this->_infDiffuseResult + this->_pointDiffuseResult +
            this->_infSpecularResult + this->_pointSpecularResult;
    return L2D::Color::colorClamp(unClamped.red, unClamped.green, unClamped.blue);
}

void L3D::LightCaster::applyTransformation(const Matrix &matrix) {
    for (L3D::Light* light : lightSources) {
        // This check has no actual meaning, as a applyTransformation() call
        // is a noop. This does make the code more explicit in what it is supposed to do.
        if (light->isPointLight() || light->isInfLight())
            light->applyTransformation(matrix);
    }
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

void L3D::Face::addToLines2DZ(L2D::Lines2DZ& lines2DZ,
                              const std::vector<Vector3D> &points,
                              const L2D::Color& color,
                              const double projectionScreenDistance) const {

    // if the face consists of only two points (a singe line),
    // then the list should not be seen as circular
    if (point_indexes.size() == 2) {
        const Vector3D &p1 = points.at(point_indexes.at(0));
        const Vector3D &p2 = points.at(point_indexes.at(1));

        lines2DZ.emplace_back(L2D::Line2DZ(projectPoint3D(p1, projectionScreenDistance), p1.z,
                                           projectPoint3D(p2, projectionScreenDistance), p2.z,
                                           color.red, color.green, color.blue));
    } else {
        // TODO if a point index pair is found in multiple faces, then we make a line for each occurrence !!! =! good
        unsigned i2 = 1;
        for (unsigned int i1 = 0; i1 < point_indexes.size(); i1++) {

            // retrieve the 3D points that need to be converted
            const Vector3D &p1 = points.at(point_indexes.at(i1));
            const Vector3D &p2 = points.at(point_indexes.at(i2));

            // project a 3D line of two 3D points onto a 2D line
            lines2DZ.emplace_back(L2D::Line2DZ(projectPoint3D(p1, projectionScreenDistance), p1.z,
                                               projectPoint3D(p2, projectionScreenDistance), p2.z,
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






L3D::Figure L3D::Figure::createLineDrawingFigure(const ini::Configuration& configuration,
                                                 const std::string& figureName) {

    L3D::Figure newFigure = L3D::Figure();

    // add all the 3D points to the L3D::Figure
    std::vector<double> point = {};
    for (unsigned int pointIndex = 0; pointIndex < (unsigned int) configuration[figureName]["nrPoints"].as_int_or_die(); pointIndex++) {

        point = configuration[figureName]["point" + std::to_string(pointIndex)];

        newFigure.points.emplace_back(Vector3D::point(point.at(0), point.at(1), point.at(2)));
    }

    // add all faces to the figure
    std::vector<int> line = {};
    for (unsigned int lineIndex = 0; lineIndex < (unsigned int) configuration[figureName]["nrLines"].as_int_or_die(); lineIndex++) {

        L3D::Face newFace;

        line = configuration[figureName]["line" + std::to_string(lineIndex)].as_int_tuple_or_die();
        newFace.point_indexes.emplace_back(line.at(0));
        newFace.point_indexes.emplace_back(line.at(1));

        newFigure.faces.emplace_back(newFace);
    }

    return newFigure;
}

L3D::Figure L3D::Figure::createCube(const bool triangulate) {
    return createBasicPlatonicBody("cube", triangulate);
}

L3D::Figure L3D::Figure::createTetrahedron() {
    return createBasicPlatonicBody("tetrahedron");
}

L3D::Figure L3D::Figure::createOctahedron() {
    return createBasicPlatonicBody("octahedron");
}

L3D::Figure L3D::Figure::createIcosahedron() {

    L3D::Figure newIcosahedron = L3D::Figure();

    ini::Configuration configuration;
    std::ifstream fin("../3D_Bodies/icosahedron.ini");
    fin >> configuration;
    fin.close();

    // generate the points for the icosahedron
    newIcosahedron.points.emplace_back(Vector3D::point(0, 0,   std::sqrt(5.0) / 2.0)); // top point
    for (unsigned int p = 2; p <= 6; p++)
        newIcosahedron.points.emplace_back(Vector3D::point(std::cos((p - 2) * 2.0 * M_PI / 5.0),
                                                           std::sin((p - 2) * 2.0 * M_PI / 5.0),
                                                           0.5));
    for (unsigned int p = 7; p <= 11; p++)
        newIcosahedron.points.emplace_back(Vector3D::point(std::cos((1.0 + ((p - 7.0) * 2.0)) * M_PI / 5.0),
                                                           std::sin((1.0 + ((p - 7.0) * 2.0)) * M_PI / 5.0),
                                                           -0.5));

    newIcosahedron.points.emplace_back(Vector3D::point(0, 0, - std::sqrt(5.0) / 2.0)); // bottom point


    // Parse the faces for the icosahedron.
    // A icosahedron is already made up of triangles,
    // sp triangulating it has no meaning.
    parseFacesPlatonicBody(configuration, newIcosahedron);

    return newIcosahedron;
}

L3D::Figure L3D::Figure::createDodecahedron(const bool triangulate) {

    L3D::Figure newDodecahedron = L3D::Figure(triangulate);

    ini::Configuration configuration;
    std::ifstream fin("../3D_Bodies/dodecahedron.ini");
    fin >> configuration;
    fin.close();

    // generate the points for the dodecahedron from a icosahedron
    const L3D::Figure icosahedron = createIcosahedron();
    Vector3D center = Vector3D::point(0, 0, 0);

    for (const L3D::Face& iFace : icosahedron.faces) {

        // Calculate the center point of the face
        for (unsigned int pointIndex : iFace.point_indexes)
            center += icosahedron.points.at(pointIndex);
        center /= 3.0;

        // The resulting center is a point of the dodecahedron
        newDodecahedron.points.emplace_back(center);

        // Reset center point variable for reuse
        center = Vector3D::point(0, 0, 0);
    }

    // parse the faces for the dodecahedron
    parseFacesPlatonicBody(configuration, newDodecahedron);

    return newDodecahedron;
}

L3D::Figure L3D::Figure::createCone(const ini::Configuration &configuration,
                                    const std::string& figureName,
                                    const bool triangulate) {

    L3D::Figure newCone = L3D::Figure(triangulate);
    const int n = configuration[figureName]["n"].as_int_or_die();
    const double height = configuration[figureName]["height"].as_double_or_die();

    if (n < 3)
        return newCone;

    // generate the points for the cone
    const double c = 2.0 * M_PI / ((double) n);
    L3D::Face newFace = L3D::Face();
    newFace.point_indexes.resize(n);    // bottom face will hold n points

    for (unsigned int p = 0; p < (unsigned int) n; p++) {
        newCone.points.emplace_back(Vector3D::point(std::cos(p * c), std::sin(p * c), 0));
        newFace.point_indexes.at(n - 1 - p) = p;    // index value decreases from left to right
    }


    // Triangulate the bottom face if needed
    newCone.addFace(newFace);

    newCone.points.emplace_back(Vector3D::point(0, 0, height)); // top point

    // generate the faces for the cone
    newFace.point_indexes.resize(3);
    const int& topIndex = n;
    for (unsigned int pointIndex = 0; pointIndex < (unsigned int) n - 1; pointIndex++) {
        newFace.point_indexes.at(0) = pointIndex;
        newFace.point_indexes.at(1) = pointIndex + 1;
        newFace.point_indexes.at(2) = topIndex;
        newCone.faces.emplace_back(newFace);
    }

    // add the last face with points [p0, pn-1, ptop]
    newFace.point_indexes.at(0) = topIndex - 1; // p0
    newFace.point_indexes.at(1) = 0;            // pn-1
    newFace.point_indexes.at(2) = topIndex;     // ptop
    newCone.faces.emplace_back(newFace);

    return newCone;
}

L3D::Figure L3D::Figure::createCylinder(const ini::Configuration &configuration,
                                        const std::string& figureName,
                                        const bool triangulate) {

    L3D::Figure newCylinder = L3D::Figure(triangulate);
    const int n = configuration[figureName]["n"].as_int_or_die();
    const double height = configuration[figureName]["height"].as_double_or_die();

    if (n < 3)
        return newCylinder;

    // generate the points for the cone
    L3D::Face newFace = L3D::Face();
    newFace.point_indexes.resize(n);
    const double c = 2.0 * M_PI / ((double) n);
    // bottom face points
    for (unsigned int p = 0; p < (unsigned int) n; p++) {
        newCylinder.points.emplace_back(Vector3D::point(std::cos(p * c), std::sin(p * c), 0));
        newFace.point_indexes.at(n - 1 - p) = p;    // index value decreases from left to right
    }
    newCylinder.addFace(newFace);

    // top face points
    for (unsigned int p = 0; p < (unsigned int) n; p++) {
        newCylinder.points.emplace_back(Vector3D::point(newCylinder.points.at(p).x, newCylinder.points.at(p).y, height));
        newFace.point_indexes.at(p) = p + n; // index value increases from left to right
    }
    newCylinder.addFace(newFace);

    // generate the faces for the cone<<
    newFace.point_indexes.resize(4);
    for (unsigned int pointIndex = 0; pointIndex < ((unsigned int) n) - 1; pointIndex++) {
        newFace.point_indexes.at(0) = pointIndex;           // pi           bottom left
        newFace.point_indexes.at(1) = pointIndex + 1;       // pi+1         bottom right
        newFace.point_indexes.at(2) = pointIndex + n + 1;   // pi+n+1       top right
        newFace.point_indexes.at(3) = pointIndex + n;       // pi+n         top left

        newCylinder.addFace(newFace);
    }

    // add the last face with points [p0, pn-1, pn+1]
    newFace.point_indexes.at(0) = n - 1;        // pi           bottom left,   last point of bottom face
    newFace.point_indexes.at(1) = 0;            // pi+1         bottom right,  first point of bottom face
    newFace.point_indexes.at(2) = n;            // pi+n+1       top right,     first point of top face
    newFace.point_indexes.at(3) = (2*n) - 1;    // pi+n         top left,      last point of top face

    newCylinder.addFace(newFace);

    return newCylinder;
}

L3D::Figure L3D::Figure::createSphere(const ini::Configuration &configuration,
                                      const std::string& figureName) {

    L3D::Figure newSphere = L3D::Figure();
    L3D::Figure icosahedron = L3D::Figure::createIcosahedron();

    const unsigned int iterationsLeft = configuration[figureName]["n"].as_int_or_die();

    // copy over the point of the icosahedron
    for (const Vector3D& point : icosahedron.points) {
        newSphere.points.emplace_back(point);
    }

    // Generate all the points of the circle
    unsigned int Ai;
    unsigned int Bi;
    unsigned int Ci;

    for (const Face& face : icosahedron.faces) {
        Ai = face.point_indexes.at(0);
        Bi = face.point_indexes.at(1);
        Ci = face.point_indexes.at(2);

        divisionTriangleFace(iterationsLeft,
                             Ai, Bi, Ci,
                             icosahedron.points.at(Ai), icosahedron.points.at(Bi), icosahedron.points.at(Ci),
                             newSphere);
    }

    for (Vector3D& point : newSphere.points) {
        point.normalise();
    }

    return newSphere;
}

L3D::Figure L3D::Figure::createBasicPlatonicBody(const std::string& type,
                                                 const bool triangulate) {

    // return an empty L3D::Figure
    if (type != "cube" && type != "tetrahedron" && type != "octahedron")
        return L3D::Figure();

    L3D::Figure newPlatonicBody = L3D::Figure(triangulate);
    ini::Configuration configuration;

    std::ifstream fin("../3D_Bodies/" + type + ".ini");
    fin >> configuration;
    fin.close();

    // add all the 3D points to the L3D::Figure
    parsePointsPlatonicBody(configuration, newPlatonicBody);

    // add all faces to the figure
    parseFacesPlatonicBody(configuration, newPlatonicBody);

    return newPlatonicBody;
}




void L3D::Figure::parsePointsPlatonicBody(const ini::Configuration& configuration, L3D::Figure& platonicBody) {

    std::vector<double> point = {};
    for (unsigned int pointIndex = 0; pointIndex < (unsigned int) configuration["Points"]["nrPoints"].as_int_or_die(); pointIndex++) {

        point = configuration["Points"]["point" + std::to_string(pointIndex)].as_double_tuple_or_die();

        platonicBody.points.emplace_back(Vector3D::point(point.at(0), point.at(1), point.at(2)));
    }
}

void L3D::Figure::triangulateFace(const std::vector<int>& pointIndexes, L3D::Figure &figure) {

    L3D::Face newFace = L3D::Face();
    newFace.point_indexes.resize(3);

    // Choose the first point in the face as the common point between
    // all generated triangular faces.
    // Then, divide the face into triangles based on that base point.
    newFace.point_indexes.at(0) = pointIndexes.at(0);
    for (unsigned int i = 1; i+1 < pointIndexes.size(); i++) {
        newFace.point_indexes.at(1) = pointIndexes.at(i);
        newFace.point_indexes.at(2) = pointIndexes.at(i+1);

        figure.faces.emplace_back(newFace);
    }
}

void L3D::Figure::parseFacesPlatonicBody(const ini::Configuration& configuration, L3D::Figure& platonicBody) {

    std::vector<int> face = {};
    for (unsigned int faceIndex = 0; faceIndex < (unsigned int) configuration["Faces"]["nrFaces"].as_int_or_die(); faceIndex++) {

        L3D::Face newFace;

        face = configuration["Faces"]["face" + std::to_string(faceIndex)].as_int_tuple_or_die();

        newFace.point_indexes = face;
        platonicBody.addFace(newFace);
    }
}

void L3D::Figure::divisionTriangleFace(const unsigned int iterationsLeft,
                                   const unsigned int Ai, const unsigned int Bi, const unsigned int Ci,
                                   const Vector3D &A, const Vector3D &B, const Vector3D &C,
                                   L3D::Figure &sphere) {

    // Create a face representing the triangle ABC
    if (!iterationsLeft) {
        L3D::Face ABC = L3D::Face();
        ABC.point_indexes.emplace_back(Ai);
        ABC.point_indexes.emplace_back(Bi);
        ABC.point_indexes.emplace_back(Ci);
        sphere.faces.emplace_back(ABC);
    }
    // Division the triangle ABC into four new triangles.
    // Then recurse on those triangles.
    else {
                                            //             C
        const Vector3D D = (A + B) / 2.0;   //
        const Vector3D E = (A + C) / 2.0;   //         E      F
        const Vector3D F = (B + C) / 2.0;   //
                                            //      A     D     B

        sphere.points.emplace_back(D);
        sphere.points.emplace_back(E);
        sphere.points.emplace_back(F);

        unsigned int pointsSize = sphere.points.size();

        const unsigned int Di = pointsSize - 3;
        const unsigned int Ei = pointsSize - 2;
        const unsigned int Fi = pointsSize - 1;

        // A structured array of the indexes of the four sub triangles
        unsigned int triangleIndexes[12] = {
                Di, Fi, Ei,
                Ai, Di, Ei,
                Di, Bi, Fi,
                Ei, Fi, Ci
        };
        // A structured array of the vertexes of the four sub triangles
        const Vector3D* triangleVertexes[12] = {
                &D, &F, &E,
                &A, &D, &E,
                &D, &B, &F,
                &E, &F, &C
        };

        for (unsigned int i = 0; i < 12; i+=3) {
            divisionTriangleFace(iterationsLeft - 1,
                                 triangleIndexes[i],  triangleIndexes[i+1],   triangleIndexes[i+2],
                                 *triangleVertexes[i], *triangleVertexes[i+1], *triangleVertexes[i+2],
                                 sphere);
        }
    }

}

L3D::Figure L3D::Figure::createTorus(const ini::Configuration &configuration,
                                     const std::string &figureName, const bool triangulate) {

    L3D::Figure newTorus = L3D::Figure(triangulate);

    // Distance from center of the hole of the torus to center of the tube of the torus.
    double R = configuration[figureName]["R"].as_double_or_die();
    // Radius of tube of the torus.
    double r = configuration[figureName]["r"].as_double_or_die();
    // The amount of vertical circles the torus is divided into.
    unsigned int n = configuration[figureName]["n"].as_int_or_die();
    // The amount of points each vertical circle is approximated with.
    // Alternately, the amount of "rectangles" between two vertical circles.
    unsigned int m = configuration[figureName]["m"].as_int_or_die();

    double cn = 2.0 * M_PI / ((double) n);
    double cm = 2.0 * M_PI / ((double) m);
    double u;   // The angle that indicates a vertical circle
    double v;   // The angle that indicates a point on a vertical circle

    // We calculate all points, one vertical circle at a time.
    // That way the points needed to form the "rectangles" between
    // neighbouring vertical circles are close together in the
    // points list and in memory.
    for (unsigned int i = 0; i < n; i++) {
        u = i * cn;
        for (unsigned int j = 0; j < m; j++) {
            v = j * cm;
            newTorus.points.emplace_back(Vector3D::point((R + r*std::cos(v))*std::cos(u),
                                                         (R + r*std::cos(v))*std::sin(u),
                                                         r * std::sin(v)));
        }
    }

    // generate the faces
    L3D::Face newFace = L3D::Face();
    newFace.point_indexes.resize(4);
    unsigned int fcircBaseInd = 0;  // The index of the first point of the first circle (vertical circle i)
    unsigned int scircBaseInd = n;  // The index of the first point of the second circle (vertical circle i+1)
    for (unsigned int circleIndex = 0; circleIndex < n; circleIndex++) {
        for (unsigned int pointIndex = 0; pointIndex < m; pointIndex++) {
            newFace.point_indexes.at(0) = fcircBaseInd + pointIndex;      // p_{ i  , j   }
            newFace.point_indexes.at(1) = scircBaseInd + pointIndex;      // p_{ i+1, j   }
            newFace.point_indexes.at(2) = scircBaseInd + ((pointIndex != m-1) ? (pointIndex + 1) : 0) ;  // p_{ i+1, j+1 }
            newFace.point_indexes.at(3) = fcircBaseInd + ((pointIndex != m-1) ? (pointIndex + 1) : 0);  // p_{ i  , j+1 }

            newTorus.addFace(newFace);
        }

        fcircBaseInd = scircBaseInd;
        if (circleIndex == m-2) // last loop: circles n-1 and 0 instead of i and i+1
            scircBaseInd = 0;
        else
            scircBaseInd += n;
    }

    return newTorus;
}


L3D::Figure::Figure(const bool triangulate) : triangulate(triangulate) {}

L3D::Figure::Figure(L2D::Color &ambientReflectivity,
                    L2D::Color &diffuseReflectivity,
                    L2D::Color &specularReflectivity,
                    double reflectionCoefficient, bool triangulate) : triangulate(triangulate) {

    this->ambientReflectivity = ambientReflectivity;
    this->diffuseReflectivity = diffuseReflectivity;
    this->specularReflectivity = specularReflectivity;
    this->reflectionCoefficient = reflectionCoefficient;
}


bool L3D::Figure::addFace(const L3D::Face &newFace) {

    if (triangulate && newFace.point_indexes.size() >= 3) {
        triangulateFace(newFace.point_indexes, *this);
        return true;
    }

    faces.emplace_back(newFace);
    return false;
}

void L3D::Figure::applyTransformation(const Matrix &M) {

    for (Vector3D& point : points) {
        point *= M;
    }
}

L2D::Lines2D L3D::Figure::toLines2D(const double projectionScreenDistance) const {

    L2D::Lines2D lines2D = {};

    for (const L3D::Face& face : faces) {
        face.addToLines2D(lines2D, points, this->ambientReflectivity, projectionScreenDistance);
    }

    return lines2D;
}

L2D::Lines2DZ L3D::Figure::toLines2DZ(const double projectionScreenDistance) const {

    L2D::Lines2DZ lines2DZ = {};

    for (const L3D::Face& face : faces) {
        face.addToLines2DZ(lines2DZ, points, this->ambientReflectivity, projectionScreenDistance);
    }

    return lines2DZ;
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




L3D::ZBuffer::ZBuffer(const unsigned int width, const unsigned int height) {
    this->width = width;
    this->height = height;

    this->resize(width);

    for (std::vector<double>& column : *this) {
        column.resize(height, +std::numeric_limits<double>::infinity());
    }
}

bool L3D::ZBuffer::shouldReplace(const unsigned int x, const unsigned int y, const double zInv) const {
    assert(x < width && y < height);
    return zInv < (*this)[x][y];
}

bool L3D::ZBuffer::replace(const unsigned int x, const unsigned int y, const double zInv) {
    if (shouldReplace(x, y, zInv)) {
        (*this)[x][y] = zInv;
        return true;
    }
    return false;
}



L3D::LSystem::State::State() : canReloadSafely(false) {
    currentLoc = Vector3D::point(0, 0, 0);
    H = Vector3D::vector(1, 0, 0);
    L = Vector3D::vector(0, 1, 0);
    U = Vector3D::vector(0, 0, 1);
}

L3D::LSystem::State::State(const double x, const double y, const double z) : canReloadSafely(false) {
    currentLoc = Vector3D::point(x, y, z);
    H = Vector3D::vector(1, 0, 0);
    L = Vector3D::vector(0, 1, 0);
    U = Vector3D::vector(0, 0, 1);
}

void L3D::LSystem::State::rotateAroundU(double angle) {
    double rAngle = toRadians(angle);
    H = H *  std::cos(rAngle)  +  L * std::sin(rAngle);
    L = H * -std::sin(rAngle)  +  L * std::cos(rAngle);
}
void L3D::LSystem::State::yaw(double angle) {
    double rAngle = toRadians(angle);
    Vector3D Hc = H;

    H = H  * std::cos(rAngle)   +  L * std::sin(rAngle);
    L = Hc * std::sin(-rAngle)  +  L * std::cos(rAngle);
}

void L3D::LSystem::State::rotateAroundL(double angle) {
    double rAngle = toRadians(angle);
    H = H *  std::cos(rAngle)  +  U * std::sin(rAngle);
    U = H * -std::sin(rAngle)  +  U * std::cos(rAngle);
}
void L3D::LSystem::State::pitch(double angle) {
    double rAngle = toRadians(angle);
    Vector3D Hc = H;

    H = H  *  std::cos(rAngle)  +  U * std::sin(rAngle);
    U = Hc *  std::sin(-rAngle) +  U * std::cos(rAngle);
}

void L3D::LSystem::State::rotateAroundH(double angle) {
    double rAngle = toRadians(angle);
    L = L * std::cos(rAngle)  -  U * std::sin(rAngle);
    U = L * std::sin(rAngle)  +  U * std::cos(rAngle);
}
void L3D::LSystem::State::roll(double angle) {
    double rAngle = toRadians(angle);
    Vector3D Lc = L;

    L = L  * std::cos(rAngle)  +  U * std::sin(-rAngle);
    U = Lc * std::sin(rAngle)  +  U * std::cos(rAngle);
}

void L3D::LSystem::State::backFlip() {
    H *= -1.0;
    L *= -1.0;
}

void L3D::LSystem::State::move() {
    currentLoc += H;
}




L3D::LSystem::LGenerator::LGenerator() : _angle(0), _justTunneled(false), _justReloaded(false) {}

L3D::Figure L3D::LSystem::LGenerator::generateFigure(const ini::Configuration &configuration,
                                                     const LParser::LSystem3D &lSystem) {

    L3D::Figure L3DFigure = L3D::Figure();

    addFaces(lSystem, L3DFigure);

    return L3DFigure;
}

void L3D::LSystem::LGenerator::addFaces(const LParser::LSystem3D &lSystem,
                                        L3D::Figure &L3DFigure) {

    _angle = lSystem.get_angle();

    _p1 = Vector3D::point(0.0, 0.0, 0.0);
    _p2 = Vector3D::point(0.0, 0.0, 0.0);

    L3DFigure.points.emplace_back(_p1);

    // The stack should be empty after each lines generation,
    // but we empty it to avoid accidental large memory consumption
    if (!_savedStates.empty())
        _savedStates = {};

    for (const char c : lSystem.get_initiator())
    {
        recurse(c, lSystem, lSystem.get_nr_iterations(), L3DFigure);
    }

}

void L3D::LSystem::LGenerator::recurse(char replacedChar,
                                       const LParser::LSystem3D &lSystem,
                                       unsigned int iterations,
                                       L3D::Figure &L3DFigure) {
    // change the angle as needed
    if (replacedChar == '+') {
        _state.yaw(_angle);
    } else if (replacedChar == '-') {
        _state.yaw(-_angle);
    } else if (replacedChar == '^') {
        _state.pitch(_angle);
    } else if (replacedChar == '&') {
        _state.pitch(-_angle);
    } else if (replacedChar == '\\') {
        _state.roll(_angle);
    } else if (replacedChar == '/') {
        _state.roll(-_angle);
    } else if (replacedChar == '|') {
        _state.backFlip();
    } else if (replacedChar == '(') {

        _state.canReloadSafely = !_justTunneled;
        // in the case of ')(', _state.savedIndex must remain unchanged
        if (!_justReloaded && _state.canReloadSafely)
            _state.savedIndex = L3DFigure.points.size()-1;

        _savedStates.push(_state);
    } else if (replacedChar == ')') {
        if (_savedStates.empty())
            std::cerr << "During the generation of a 3D LSystem, a ')' was encountered."
                      << "There was, however, no saved stated to reload to." << std::endl;

        _state = _savedStates.top();
        _savedStates.pop();

        // only allow the use of savedIndex if it is safe to do so
        _justReloaded = _state.canReloadSafely;
        _justTunneled = false;
    }
        // if recursion necessary, do it
    else if (iterations > 0)
    {
        for (const char c : lSystem.get_replacement(replacedChar))
        {
            recurse(c, lSystem, iterations-1, L3DFigure);
        }
    }
        // no more recursion needed, so start adding lines
    else if (lSystem.draw(replacedChar)) {

        // starting point of the new line
        _p1 = _state.currentLoc;

        // adjust the current coordinates
        _state.move();

        // end point of the new line
        _p2 = _state.currentLoc;

        // add a new L3D::Face to the figure
        if (_justTunneled)
            L3DFigure.points.emplace_back(_p1);     // TODO duplicates !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        L3DFigure.points.emplace_back(_p2);

        L3D::Face newLine = L3D::Face();
        newLine.point_indexes.emplace_back(_justReloaded ? _state.savedIndex : L3DFigure.points.size()-2);
        newLine.point_indexes.emplace_back(L3DFigure.points.size()-1);

        L3DFigure.faces.emplace_back(newLine);

        _justTunneled = false;
        _justReloaded = false;

    } else {

        // adjust the current coordinates
        _state.move();

        _justTunneled = true;
        _justReloaded = false;
    }


}

