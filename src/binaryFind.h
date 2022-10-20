/**
 * @file binaryFind.cpp
 *
 * @brief Fundamental Logic Class: Simple set operations (Is in set? etc.)
 *
 * @ingroup Fundamental
 *
 *
 * @author Chris Prior
 * Contact: christopher.prior@durham.ac.uk
 *
 */
#ifndef BINARY_FIND_H
#define BINARY_FIND_H

#include <algorithm>
#include "point.h"

class binaryFind{
	public:
	binaryFind();
	/**
	 * Checks if input value lies within a given interval. Has no return value.
	 * 
	 * @param lowerSec Minimum of interval.
	 * @param upperSec Maximum of interval.
	 * @param zval Input Value.
	 */
	void checkInSec(double lowerSec, double upperSec,double& zval);
	/**
	 * Checks if input value lies within a given interval.
	 * Do we need both?
	 * 
	 * @param lowerSec Minimum of interval.
	 * @param upperSec Maximum of interval.
	 * @param zv Input Value.
	 * @return true 
	 * @return false 
	 */
	bool checkInSecInitial(double lowerSec, double upperSec,double& zv);
	/**
	 * Gets the pair of points whose z-value bounds the input value.
	 * 
	 * @param pointList List of points.
	 * @param lowerIndex Index of first point.
	 * @param upperIndex Index of last point.
	 * @param zv Input value.
	 * @return std::pair<int,int> 
	 */
	std::pair<int,int> getContainingPair(std::vector<point>& pointList,int& lowerIndex,int& upperIndex,double& zv);
	private:
	bool isIn;	
	int lowerIndex,upperIndex;
	std::pair<int,int> boundingIndicies;
	double zval;
};

#endif
