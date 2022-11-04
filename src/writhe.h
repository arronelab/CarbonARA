#ifndef WRITHE_H
#define WRITHE_H

#include "point.h"
#include <complex>

class writhe{
public:
    writhe();
    /**
     * Computes the vector between the ith and jth point of pointList.
     * 
     * @param pointList List of coordinates
     * @param size Unsure, ith and jth coordinate are taken modulo this paramater. I guess ensures we stay in the pointlist?
     * @param i Index of start point
     * @param j Index of end point
     * @return Vector
     */
    point dbold(std::vector<point>& pointList,int size,int i,int j);
    /**
     * Computes the unit vector between the ith and (i-1)th point of pointList.
     * 
     * @param pointList List of coordinates
     * @param size Unsure, ith and jth coordinate are taken modulo this paramater. I guess ensures we stay in the pointlist?
     * @param i Index of end point.
     * @return Vector 
     */
    point te(std::vector<point>& pointList,int size,int i);
    /**
     * Computes inverse cosine of input. This is useful in giving the angle the between vectors, which is a component of the writhe calculation.
     * 
     * @param temp Input value.
     * @return Angle.
     */
    double acosC(double temp);
    double mu(std::vector<point>& pointList,int size,int i,int k,int m, int j);
    /**
     * Computes Omega_ij as per https://en.wikipedia.org/wiki/Writhe subsection: Numerically approximating the Gauss integral for writhe of a curve in space
     * 
     * @param pointList List of coordinates
     * @param size Unsure, ith and jth coordinate are taken modulo this paramater. I guess ensures we stay in the pointlist?
     * @param i Start index
     * @param j End index
     * @return double 
     */
    double wij(std::vector<point>& pointList,int size,int i,int j);
    /**
     * Computes the Average Crossing Number (ACN) for an input curve via numerical approximation of Gauss Linking Integral.
     * 
     * @param pointList Coordinates of discrete curve.
     * @return double 
     */
    double DIAbs(std::vector<point>& pointList);
    /**
     * Takes an input curve (should be a list of lists, where each sublist is a subsection of the curve, think secondary structural elements)
     * Smooths the curve so that each edge represents a subsection of the input curve,
     * Computes the ACN of this smoothed curve.
     * 
     * @param pointListIn Input curve, should be a list of lists.
     * @return double 
     */
    double DIDownSampleAbs(std::vector<std::vector<point> >& pointListIn);
    /**
     * Takes an input curve (should be a list of lists, where each sublist is a subsection of the curve, think secondary structural elements)
     * Smooths the curve so that each edge represents a subsection of the input curve,
     * Computes the ACN FingerPrint of this smoothed curve.
     * 
     * @param pointListIn Input curve, should be a list of lists.
     * @return std::vector<std::pair<std::pair<int,int>,double> > 
     */
    std::vector<std::pair<std::pair<int,int>,double> > DIDownSampleAbsFP(std::vector<std::vector<point> >& pointListIn);
};

#endif