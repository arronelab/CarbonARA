#include "writhe.h"

const double PI  =3.141592653589793238463;

point writhe::dbold(std::vector<point>& pointList,int size,int i,int j){
  i = i%(size-1);
  j = j%(size-1);
  return pointList[i+1].dif(pointList[j+1]);
}

point writhe::te(std::vector<point>& pointList,int size,int i){
  point temp = dbold(pointList,size,i,i-1);
  temp.normalise();
  return temp;
}

double writhe::acosC(double temp){
  if(temp < -1.0){
    return PI;
  }else if(temp > 1.0){
    return 0;
  }else{
    return acos(temp);
  }
}

double writhe::mu(std::vector<point>& pointList,int size,int i,int k,int m, int j){
  point ti= te(pointList,size,i);
  point tj=te(pointList,size,j);
  point dkm = dbold(pointList,size,k,m);
  point un = ti.cross(dkm);
  point deux = dkm.cross(tj);
  if( !un.isNonzero() || !deux.isNonzero()){
    return 0;
  }else{
    un.normalise();
    deux.normalise();
    double temp = un.dotprod(deux);
    point cp = ti.cross(tj);
    double signProd = dkm.dotprod(cp); 
    if(std::abs(signProd)<0.000000001){
      return 0.0;
    }else{
      if(signProd>= 0){
	return acosC(temp);
      }else{
	return (-1.0)*acosC(temp);
      }
    }
   }
}

double writhe::wij(std::vector<point>& pointList,int size,int i,int j){
  double temp =0;
  if(j>i){
    double t1 = mu(pointList,size,i,i-1,j-1,j);
    double t2 = mu(pointList,size,i,i,j-1,j);
    double t3 = mu(pointList,size,i,i-1,j,j);
    double t4 = mu(pointList,size,i,i,j,j);
    temp = t1-t2-t3+t4;
  }
  return temp;
}

double writhe::DIAbs(std::vector<point>& pointList){
  int listSize = pointList.size();
  double sigsum=0.0;
  for(int i=0;i<listSize-1;i++){
    for(int j=i+1;j< listSize-1;j++){
      sigsum = sigsum +std::abs(wij(pointList,listSize,i,j));
     }
  }
  return sigsum/(2*PI);
}

std::vector<std::pair<std::pair<int,int>,double> > writhe::DIDownSampleAbs(std::vector<std::vector<point> >& pointListIn){
  // first downsample
  std::vector<point> mol;
  for(int i=0;i<pointListIn.size();i++){
    mol.push_back(pointListIn[i][0]);
  }
  std::vector<point> lastSec = pointListIn[pointListIn.size()-1];
  mol.push_back(lastSec[lastSec.size()-1]);
  std::vector<std::pair<std::pair<int,int>,double> > fingerPrintList;
  for(int l=5;l<mol.size();l++){
     for(int k=0;k<mol.size()-5;k++){
       if(k+l<mol.size()){
	 // grab subsection of curve;
	 std::vector<point>::const_iterator first = mol.begin()+k;
         std::vector<point>::const_iterator last = mol.begin()+k+l;
	 std::vector<point> subVec(first, last);
	 double wrval = DIAbs(subVec);
         std::pair<int,int> pr;
         pr.first=k;pr.second=k+l;     
         std::pair<std::pair<int,int>,double> fppoint;
         fppoint.first=pr;fppoint.second=wrval;
         fingerPrintList.push_back(fppoint);
       }
    }   
  }
  return fingerPrintList;
}