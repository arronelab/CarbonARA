/**
 * @file point.cpp
 *
 * @brief Fundamental Vector Class: Standard vector operations, e.g. scalar product, cross product, euclidean norm.
 *
 * @ingroup Fundamental
 *
 *
 * @author Chris Prior
 * Contact: christopher.prior@durham.ac.uk
 *
 */
#ifndef PNT_H
#define PNT_H
#include<iostream>
#include<fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <cmath>
#include <vector>


class point{
            public:
            /**
             * Constructs a vector with given coordinates
             * 
             * @param x x-coordinate
             * @param y y-coordinate
             * @param z z-coordinate
             */
            point(double x,double y,double z);
            /**
             # @overload Constructs an empty vector
             */
            point();
            /**
             # @overload Constructs a vector with coordinates given as string
             */
            point(std::string& triplet);
            /**
             * Get x-coordinate
             * 
             * @return double  
             */
            double getX();
            /**
             * Get y-coordinate
             * 
             * @return double  
             */
            double getY();
            /**
             * Get z-coordinate
             * 
             * @return double  
             */
            double getZ();
            /**
             * Change x-coordinate
             * 
             * @param xv New x-coordinate
             */
            void setX(double xv);
            /**
             * Change x-coordinate
             * 
             * @param yv New y-coordinate
             */
            void setY(double yv);
            /**
             * Change z-coordinate
             * 
             * @param zv New z-coordinate
             */
            void setZ(double zv);
            /**
             * Get length of vector via sqrt(x^2+y^2+z^2)
             * 
             * @return double 
             */
            double length();
            /**
             * Adds input vector to current vector (component wise)
             * 
             * @param b Vector to be added
             * @return point
             */
            point sum(point& b);
            /**
             * Subtracts input vector from current vector (component wise)
             * 
             * @param b Vector to be subtracted
             * @return point
             */
            point dif(point& b);
            /**
             * Computes cross product of input vector and current vector
             * 
             * @param p2 Vector to be used in cross product
             * @return point
             */
            point cross(point& p2);
            /**
             * Computes scalar triple product of three vectors p1, p2, p3
             * 
             * @param p1 Input vector
             * @param p2 Input vector
             * @param p3 Input vector
             * @return double 
             */
            double scalarTriple(point& p1,point& p2,point& p3);
            /**
             * Checks if input vector is the same as current vector
             * 
             * @param p2 Input vector
             * @return true
             * @return false 
             */
            bool checkEqual(point& p2);
            /**
             * Gives the scalar multiple of current vector
             * 
             * @param a Scalar
             */
            void scalarMult(double a);
            /**
             * Checks if the pairwise distance between all components of input vector 
             * and current vector is less than 0.00000001 (i.e. x1-x0<0.00000001 and y1-y0<0.00000001 and z1-z0<0.00000001)
             * 
             * @param p Input vector
             * @return double 
             */
            double pairDist(point &p);
            /**
             * Normalises current vector
             */
            void normalise();
            /**
             * Rescales current vector such that z-component is 1
             * 
             */
            void znormalise();
            /**
             * Prints current point as X Y Z
             * 
             */
            void printPoint();
            /**
             * Computes dot product of input vector and current vector.
             * 
             * @param p2 Input vector
             * @return double 
             */
            double dotprod(point& p2);
            /**
             * Checks is all components of current vector are > 0.00000001
             * 
             * @return true 
             * @return false 
             */
            bool isNonzero();
            /**
             * Adds input vector to current vector (component wise)
             * This feels the same as sum, do we need both?
             * 
             * @param p Vector to be added
             * @return point
             */
            point operator+(point p);
            /**
             * Substracts input vector from current vector (component wise)
             * This feels the same as dif, do we need both?
             * 
             * @param p Vector to be subtracted
             * @return point
             */
            point operator-(point p);
            /**
             * Gives scalar multiple of current vector
             * Feels same as scalarMult, do we need both?
             * 
             * @param d Scalar
             * @return point 
             */
            point operator*(double d);
            /**
             * Gives scalar division of current vector (i.e. (x,y,z)-->(x/d,y/d,z/d) )
             * 
             * @param d Scalar
             * @return point 
             */
            point operator/(double d);
            /**
             * Computes Euclidean distance between input vector and current vector
             * 
             * @param p2 Input vector
             * @return double 
             */
            double eDist(point &p2);
      private:
            double X,Y,Z,norm;
};

#endif
