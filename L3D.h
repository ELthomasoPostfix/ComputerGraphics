//
// Created by Thomas Gueutal on 14/03/2021.
//

#ifndef ENGINE_L3D_H
#define ENGINE_L3D_H

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
             *      Convert the face to a list of 2D lines,
             *      which can be drawn onto an img::EasyImage.
             *      The points of the L3D::Figure should already
             *      be in eyePoint coordinates.
             *      Every line also stores the z-coordinates of
             *      the to project points for later use in z-buffering.
             *
             * @param lines2DZ:
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
            void addToLines2DZ(L2D::Lines2DZ& lines2DZ,
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
        /*
         * PUBLIC STATIC METHODS
         */
        public:

            /*
             * Description:
             *      Create and return a L3D::Figure based on the lineDrawing format.
             *
             * @param color:
             *      The color of each of the lines.
             * @param configuration:
             *      The configuration from which to parse the lines/faces.
             * @param figureName:
             *      The name by which to address the figure in the configuration.
             */
            static L3D::Figure createLineDrawingFigure(const L2D::Color& color,
                                                       const ini::Configuration& configuration,
                                                       const std::string& figureName);

            /*
             * Description:
             *      Create and return a L3D::Figure cube which is centered around (0, 0, 0).
             *      Each edge is of length 1.
             *      The resultant faces enumerate the involved points in counter clockwise fashion,
             *      when looking at a face from outside the cube.
             *      !!! Makes use of the input file /engine/3D_Bodies/cube.ini !!!
             *
             * @param color:
             *      The color of the lines of the cube.
             */
            static L3D::Figure createCube(const L2D::Color& color);

            /*
             * Description:
             *      Create and return a L3D::Figure tetrahedron which is centered around (0, 0, 0).
             *      The resultant faces enumerate the involved points in counter clockwise fashion,
             *      when looking at a face from outside the tetrahedron.
             *      !!! Makes use of the input file /engine/3D_Bodies/tetrahedron.ini !!!
             *
             * @param color:
             *      The color of the lines of the tetrahedron.
             */
            static L3D::Figure createTetrahedron(const L2D::Color& color);

            /*
             * Description:
             *      Create and return a L3D::Figure octahedron which is centered around (0, 0, 0).
             *      The resultant faces enumerate the involved points in counter clockwise fashion,
             *      when looking at a face from outside the octahedron.
             *      !!! Makes use of the input file /engine/3D_Bodies/octahedron.ini !!!
             *
             * @param color:
             *      The color of the lines of the octahedron.
            */
            static L3D::Figure createOctahedron(const L2D::Color& color);

            /*
             * Description:
             *      Create and return a L3D::Figure icosahedron which is centered around (0, 0, 0).
             *      The resultant faces enumerate the involved points in counter clockwise fashion,
             *      when looking at a face from outside the icosahedron.
             *      !!! Makes use of the input file /engine/3D_Bodies/icosahedron.ini !!!
             *
             * @param color:
             *      The color of the lines of the icosahedron.
            */
            static L3D::Figure createIcosahedron(const L2D::Color& color);
            
            /*
             * Description:
             *      Create and return a L3D::Figure dodecahedron which is centered around (0, 0, 0).
             *      The resultant faces enumerate the involved points in counter clockwise fashion,
             *      when looking at a face from outside the dodecahedron.
             *      !!! Makes use of the input file /engine/3D_Bodies/dodecahedron.ini !!!
             *
             * @param color:
             *      The color of the lines of the dodecahedron.
            */
            static L3D::Figure createDodecahedron(const L2D::Color& color);

            /*
             * Description:
             *      A function that parses a ini::Configuration file and returns
             *      the matching basic platonic body in the form of a L3D::Figure.
             *      The ini::configuration files are located at /engine/3D_Bodies/.
             *      The supported platonic bodies are 'cube', 'tetrahedron' and 'octahedron'.
             *
             * @param color:
             *      The color of each line of the basic platonic body.
             * @param type:
             *      The type of platonic body to create.
             *      The three supported types are 'cube', 'tetrahedron' and 'octahedron'.
             */
            static L3D::Figure createBasicPlatonicBody(const L2D::Color& color, const std::string& type);

            /*
             * Description:
             *      Create and return a L3D::Figure that approximates a cone.
             *      We approximate the mantle of the cone by dividing it into n triangles.
             *      n Equidistant points will be put onto the circumference
             *      of the bottom face, such that every pair of points, together with the
             *      top of the cone forms a triangle.
             *
             * @param color:
             *      The color of each of the lines.
             * @param configuration:
             *      The configuration from which to retrieve n and the height of the cone.
             * @param figureName:
             *      The name by which to address the figure in the configuration.
            */
            static L3D::Figure createCone(const L2D::Color &color,
                                              const ini::Configuration &configuration,
                                              const std::string& figureName);

            /*
             * Description:
             *      Create and return a L3D::Figure that approximates a cylinder.
             *      We approximate the mantle of the cylinder by dividing it into n rectangles.
             *      n Equidistant points will be put onto the circumference
             *      of the bottom and top faces, such that every two subsequent pairs
             *      of points form a rectangle.
             *
             * @param color:
             *      The color of each of the lines.
             * @param configuration:
             *      The configuration from which to retrieve n and the height of the cylinder.
             * @param figureName:
             *      The name by which to address the figure in the configuration.
            */
            static L3D::Figure createCylinder(const L2D::Color &color,
                                                  const ini::Configuration &configuration,
                                                  const std::string& figureName);


            /*
             * Description:
             *      Create and return a L3D::Figure that approximates a sphere.
             *      We approximate the sphere starting from a icosahedron. We
             *      then division each of the faces into smaller triangles n times.
             *      We finish by moving each resultant point to be at distance 1 from
             *      The center of the sphere.
             *
             * TODO This function creates many duplicate (Vector3D) points, the middle of
             *      the 4 triangles duplicates the points D, E and F. The 3 triangles surrounding
             *      each starting triangle/face duplicate the points generated on the edges of that
             *      starting triangle. Fix this inefficiency???
             *
             * @param color:
             *      The color of each of the lines.
             * @param configuration:
             *      The configuration from which to retrieve n (the nr of times to division each icosahedron face).
             * @param figureName:
             *      The name by which to address the figure in the configuration.
            */
            static L3D::Figure createSphere(const L2D::Color& color,
                                            const ini::Configuration& configuration,
                                            const std::string& figureName);


            /*
             * Description:
             *      Create and return a L3D::Figure that approximates a Torus.
             *      We approximate the torus by subdividing the torus into n
             *      vertical circles and then we division each circle into m
             *      points. We then create the mantle of the torus by connecting
             *      points of the circles i and j as (p_i, p_j, p_j+1, p_i+1).
             *
             * @param color:
             *      The color of each of the lines.
             * @param configuration:
             *      The configuration from which to retrieve n (nr of vertical circles)
             *      and m (nr of points on each vertical circle).
             * @param figureName:
             *      The name by which to address the figure in the configuration.
            */
            static L3D::Figure createTorus(const L2D::Color& color,
                                           const ini::Configuration& configuration,
                                           const std::string& figureName);


        /*
         * PRIVATE STATIC METHODS
         */
        private:

            /*
             * Description:
             *      Parse the points in the [Points] section of configuration and
             *      add them to the points member of platonicBody.
             *
             * @param configuration:
             *      The configuration which to parse the points in.
             * @param platonicBody:
             *      The body to which to add the parsed points.
             *
             * @return:
             *      If the type is supported, then the basic platonic body is returned.
             *      If the type is not supported, ta L3D::Figure without points or faces is returned.
             */
            static void parsePointsPlatonicBody(const ini::Configuration& configuration, L3D::Figure& platonicBody);

            /*
            * Description:
            *      Parse the faces in the [Faces] section of configuration and
            *      add them to the faces member of platonicBody.
            *
            * @param configuration:
            *      The configuration which to parse the faces in.
            * @param platonicBody:
            *      The body to which to add the parsed faces.
            */
            static void parseFacesPlatonicBody(const ini::Configuration& configuration, L3D::Figure& platonicBody);

            /*
             * Description:
             *      If iterationsLeft > 0, then we divide the triangle ABC into 4 sub triangles.
             *      We do this by calculating the middle of each edge of ABC, connecting those
             *      middle points and then recognizing each resulting triangle.
             *      When iterationsLeft == 0, then we create a L3D::Face out of ABC
             *      and add it to sphere.
             *
             * @param iterationsLeft:
             *      The number of times we need to division ABC into sub triangles and then recurse
             *      on the resulting sub triangles.
             * @param A:
             *      Corner vertex A of ABC.
             * @param B:
             *      Corner vertex B of ABC.
             * @param C:
             *      Corner vertex C of ABC.
             * @param sphere:
             *      The sphere to which to add all the resulting faces of
             *      the division.
             */
            static void divisionTriangleFace(unsigned int iterationsLeft,
                                         const unsigned int Ai, const unsigned int Bi, const unsigned int Ci,
                                         const Vector3D& A, const Vector3D& B, const Vector3D& C,
                                         L3D::Figure& sphere);

        /*
         * NON STATIC METHODS
         */
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
             *      Convert the L3D::Figure to a list of projected 2D lines.
             *      Additionally store the z-coordinates of the end points of
             *      the projected line for later use in z-buffering.
             *      We project the 3D lines onto a projection screen/plane.
             *      The screen/plane is located at a distance of projectionScreenDistance
             *      from the eyePoint and the connecting line between the
             *      eyePoint and the origin from before the eyePoint transformation
             *      was applied is perpendicular to that plane.
             *
             * @param projectionScreenDistance:
             *      The distance to the eyePoint the projection screen is located.
             */
            L2D::Lines2DZ toLines2DZ(double projectionScreenDistance) const;

            /*
             * Description:
             *      Convert the faces that make up the face into strings
             *      and list them.
             */
            std::string toString() const;

            std::ostream & operator <<(std::ostream& output_stream) const;

        public:

            std::vector<Vector3D> points;
            std::vector<Face> faces;
            L2D::Color color;
    };


    /*
     * Description:
     *      A buffer of doubles. Intended to store the smallest/lowest negative 1/z value seen yet.
     *      For that reason each element in the buffer will be initialised to +infinity.
     *
     * @member width:
     *      The width of the L3D::ZBuffer. This also denotes the
     *      amount of columns in the L3D::ZBuffer.
     * @member height:
     *      The height of the L3D::ZBuffer. This also denotes the
     *      length of each column vector.
     */
    class ZBuffer : public std::vector<std::vector<double>> {

        public:

            /*
             * Description:
             *      Construct the L3D::ZBuffer with width columns, each
             *      of which is of length height.
             *
             * @param width:
             *      The width (amount of columns) of the ZBuffer.
             *      Mirrors the width of the img::EasyImage for
             *      which the buffer is used.
             *
             * @param height:
             *      The height (length of each column) of the ZBuffer.
             *      Mirrors the width of the img::EasyImage for which
             *      the buffer is used.
             */
            ZBuffer(unsigned int width, unsigned int height);

            /*
             * Description:
             *      Determine whether zInv <  L3D::ZBuffer[x][y].
             *
             * @param x:
             *      Which column of the L3D::ZBuffer to address.
             * @param y:
             *      Which element of column x of the L3D::ZBuffer to address.
             * @param zInv:
             *      The value for which to test  L3D::ZBuffer[x][y] < zInv.
             */
            inline bool shouldReplace(unsigned int x, unsigned int y, double zInv) const;

            /*
             * Description:
             *      Replace L3D::ZBuffer[x][y] if zInv < L3D::ZBuffer[x][y].
             *      Return whether or not the replacement succeeded.
             *
             * @param x:
             *      Which column of the L3D::ZBuffer to address.
             * @param y:
             *      Which element of column x of the L3D::ZBuffer to address.
             * @param zInv:
             *      The value for which to test L3D::ZBuffer[x][y] < zInv.
             */
            bool replace(unsigned int x, unsigned int y, double zInv);

        private:

            unsigned int width;
            unsigned int height;
    };


    namespace LSystem {

        class LGenerator;

        /*
         * Description:
         *       This struct represents the mutable, base components of the 3D LGenerator.
         *       They support the generation process.
         *
         * @member currentLoc:
         *      The current x-coordinate, y-coordinate and z-coordinate in 3D space.
         * @member H:
         *      The current direction in which to draw.
         *      Initialised at (1, 0, 0).
         * @member L:
         *      The direction we consider to be "left".
         *      L is perpendicular to H.
         *      Initialised at (0, 1, 0).
         * @member U:
         *      The direction we consider to be "up".
         *      U is perpendicular to H, and is determined
         *      once we have decided H and L.
         *      Initialised at (0, 0, 1).
         */
        struct State {
            public:

                State();

                State(double x, double y, double z);

                /*
                 * Description:
                 *      A yaw rotation, which corresponds to a rotation
                 *      around the U axis.
                 *
                 * @param angle:
                 *      The angle over which to rotate in degrees.
                 */
                void rotateAroundU(double angle);
                void yaw(double angle);

                /*
                 * Description:
                 *      A pitch rotation, which corresponds to a rotation
                 *      around the L axis.
                 *
                 * @param angle:
                 *      The angle over which to rotate in degrees.
                 */
                void rotateAroundL(double angle);
                void pitch(double angle);

                /*
                 * Description:
                 *      A roll rotation, which corresponds to a rotation
                 *      around the H axis.
                 *
                 * @param angle:
                 *      The angle over which to rotate in degrees.
                 */
                void rotateAroundH(double angle);
                void roll(double angle);

                /*
                 * Description:
                 *      Perform a pitch rotation over and angle of 180°.
                 *      More performant than pitch(180).
                 */
                void backFlip();

                /*
                 * Description:
                 *      Adjust the current coordinates according to the current
                 *      direction as described by H, L and U. In other words,
                 *      move one unit.
                 */
                void move();

            public:
                Vector3D currentLoc;

            private:
                Vector3D H;
                Vector3D L;
                Vector3D U;

            protected:
                friend LGenerator;
                unsigned int savedIndex;
                bool canReloadSafely;
        };



        /*
         *
         * Description:
         *      The LGenerator class takes in a ini::Configuration and a LParser::LSystem3D.
         *      The generator creates and returns a L3D::Figure.
         *      The configuration is used to construct the returned L3D::Figure.
         *      The lsystem is used to generate the list of L3D::Faces to be added to the
         *      resulting L3D::Figure.
         *
         * @member _state:
         *      The current location and direction, represented as the state of the LGenerator.
         *      The state contains the all mutable, base attributes that support the generation process.
         * @member angle:
         *      The angle based on which to modify the current direction any time it is needed.
         *      !! Stored in degrees !!
         * @member _p1:
         *      The start point of any to generate line/L3D::Face.
         *      _p1 is (re)assigned just before adjusting the current coordinates.
         * @member _p2:
         *      The end point of any to generate line/L3D::Face.
         *      _p2 is (re)assigned just after adjusting the current coordinates.
         * @member _savedStates:
         *      It is the set of states saved during the generation of the lines/L3D::Faces.
         *      They are only valid within the generation they were saved in.
         *      If a '(' is encountered, then _state will be saved on the stack.
         *      When the corresponding ')' is encountered, then the current state
         *      will be set to the last saved state.
         */

        class LGenerator {
            private:

                /*
                 * Description:
                 *      Use the initiator of the provided lSystem to generate the lines
                 *      to be drawn, in the form of L3D::Faces.
                 *
                 * @param lSystem:
                 *      The lSystem based on which to generate the lines.
                 * @param lines:
                 *      The list to which to add the lines.
                 */
                 void addFaces(const LParser::LSystem3D& lSystem,
                               L3D::Figure& L3DFigure);


                /*
                 * Description:
                 *      Based on the replacedChar (not '+', '-', '^', '&', '\', '/', '(', ')'),
                 *      if the nr of iterations left is > 0, then recurse onto the
                 *      replacement rule of replacedChar.
                 *      Else, adjust the current direction ('+', '-', '^', '&', '\', '/'),
                 *      save or load the current state ('(', ')') or possibly generate a lineLL3D::Face.
                 *
                 * @param replacedChar:
                 *      The character on which we are recursing.
                 * @param lSystem:
                 *      The lSystem we retrieve the replacement rules from.
                 * @param iterations
                 *      The nr of recursion iterations left.
                 * @param L3DFigure:
                 *      The L3D::Figure to which to add the generated lines/L3D::Faces.
                */
                void recurse(char replacedChar,
                             const LParser::LSystem3D& lSystem,
                             unsigned int iterations,
                             L3D::Figure& L3DFigure);

        public:

                LGenerator();


                /*
                 * Description:
                 *      Generate L3D::Figure based on the ini::Configuration.
                 *
                 * @param configuration:
                 *      The requirements for the to generate L3D::Figure.
                 * @param lSystem:
                 *      The lSystem based on which to generate the lines/L3D::Faces to draw
                 *      of the to generate the L3D::Figure.
                 */
                L3D::Figure generateFigure(const ini::Configuration &configuration,
                                           L2D::Color color,
                                           const LParser::LSystem3D& lSystem);


        private:

                L3D::LSystem::State _state;

                double _angle{};

                Vector3D _p1;
                Vector3D _p2;

                std::stack<L3D::LSystem::State> _savedStates;

                bool _justTunneled;
                bool _justReloaded;
        };


    }


};


#endif //ENGINE_L3D_H
