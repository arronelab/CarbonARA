/**
 * @file optimizationPoint.cpp
 *
 * @brief Fundamental Vector Class: Processes list of points with their curvatures, torsions and lengths.
 * Unused ?? "everything's a helix" POV
 *
 * @ingroup Fundamental
 *
 *
 * @author Chris Prior
 * Contact: christopher.prior@durham.ac.uk
 *
 */
#ifndef OPTPOINT_H
#define OPTPOINT_H

#include "point.h"
#include <vector>

class optimizationPoint{
public:
/**
 * Initialises a new optimization Point object.
 * 
 * @param pointList List of points
 * @param kderivList List of curvature derivatives (?)
 * @param tderivList List of torsion derivatives (?)
 * @param LderivList List of length derivatives (?)
 */
optimizationPoint(std::vector<point> &pointList,std::vector<std::vector<point> > &kderivList,std::vector<std::vector<point> > &tderivList,std::vector<std::vector<point> > &LderivList);
/**
 * Gets the list of points from the optimization point object.
 * 
 * @return std::vector<point> 
 */
std::vector<point> getPoints();
/**
 * Gets the list of curvature derivatives from the optimization point object.
 * 
 * @return std::vector<std::vector<point> > 
 */
std::vector<std::vector<point> > getKderivs();
/**
 * Gets the list of torsion derivatives from the optimization point object.
 * 
 * @return std::vector<std::vector<point> > 
 */
std::vector<std::vector<point> > getTderivs();
/**
 * Gets the list of length derivatives from the optimization point object.
 * 
 * @return std::vector<std::vector<point> > 
 */
std::vector<std::vector<point> > getLderivs();
private:
std::vector<point> points;
std::vector<std::vector<point> > kderivs;
std::vector<std::vector<point> > tderivs;
std::vector<std::vector<point> > Lderivs;
};

#endif
