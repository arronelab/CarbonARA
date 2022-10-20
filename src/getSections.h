/**
 * @file getSections.cpp
 *
 * @brief Fundamental Logic Class: Finds turning points, necessary for localWrithe
 *
 * @ingroup Fundamental
 *
 *
 * @author Chris Prior
 * Contact: christopher.prior@durham.ac.uk
 *
 */
#ifndef GET_SECTIONS_H
#define GET_SECTIONS_H

#include "point.h"

class getSections{
      public:
      /**
       * Initialises list
       * 
       * @param lst list
       */
    	 getSections(std::vector<point>& lst);
      /**
       * Checks if there is a turning point between n1 and n2, if so, adds their positions to initialised list.
       * 
       * @param n1 Point 1
       * @param n1pos Index of point 1
       * @param n2 Point 2
       * @param n2pos Index of point 2
       * @param n3 Unused
       * @param n3pos Unused
       */
    	 void checkTurn(point& n1,int n1pos,point& n2,int n2pos,point& n3,int n3pos);
        /**
         * Finds all the turning points of a given list.
         * 
         * @return std::vector<int> List of indices of turning points of input list.
         */
         std::vector<int> getListTurnPts();
      private: 
         std::vector<point> list;
         double dif1,dif2,prod;
         std::vector<int> listSections /*think this needs to be intialised as empty*/;
      };
      
#endif
