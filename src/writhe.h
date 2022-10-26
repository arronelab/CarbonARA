#ifndef WRITHE_H
#define WRITHE_H

#include "point.h"
#include <complex>

class writhe{
public:
    point dbold(std::vector<point>& pointList,int size,int i,int j);
    point te(std::vector<point>& pointList,int size,int i);
    double acosC(double temp);
    double mu(std::vector<point>& pointList,int size,int i,int k,int m, int j);
    double wij(std::vector<point>& pointList,int size,int i,int j);
    double DIAbs(std::vector<point>& pointList);
    std::vector<std::pair<std::pair<int,int>,double> > DIDownSampleAbs(std::vector<std::vector<point> >& pointListIn);
};

#endif