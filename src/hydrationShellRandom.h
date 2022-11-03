#ifndef HYD_SH_M
#define HYD_SH_M

#include "ktlMoleculeRandom.h"
#include <algorithm>  
#include <functional> 


class hydrationShellMinimal{
public:
  hydrationShellMinimal(){};
  hydrationShellMinimal(ktlMolecule &molIn,double RinInp,double RoutInp,double &RhyIn,int ntrivsIn,double helixRatioIn,int solventsPerLinkIn,double &mutualDistCutOffIn,double &rmin,double &rmax,double &lmin);
  void getPointAndMidlengthMulti(int i,int &hIndex);
  void getPointAndMidlengthStraight(int &sec,int &part,int &hIndex,std::string soe);
  void tubeParamList();
  void getAllHelices();
  point getDirec(int i);
  point getHydTan(int i);
  point getHydNorm(int i);
  point getHydBinorm(int i);
  double hydTubeLength(int i);
  point hydTubeCentre(int i);
  point getTangent(int index,int subindex);
  point getNormal(int index,int subindex);
  point getBinormal(int index,int subindex);
  point getCentrePoint(int index,int subIndex);
  int getNScat();
  void resetMolecule();
  int getMolSize();
  int getNoSections();
  int getNoKvals();
  double getMaxDist();
  void makeInitialSegData(point &cp,point &T,point &N1,double &tm,int index,int &nseg);
  void constructInitialState();
  void generateHydrationLayer();
  std::vector<std::vector<point> > getHydrationLayer(int i);
  void solventMoleculeDistances(std::vector<double> &molSolDistances,std::vector<double> &solSolDistances);
   std::vector<std::vector<point> > returnFlatSolList();
   void writeHydrationShellToFile(const char* filename);
 private:
  ktlMolecule mol;
  std::vector<point> direcList;
  std::vector<point> frameTan;
  std::vector<point> frameNorm;
  std::vector<point> frameBinorm;
  std::vector<point> midPointList;
  std::vector<std::string> nameLst;
  std::vector<int> lengthSec;
  std::vector<double> halfLengthList,cstepList;
  int molsize;
  double Pi,Rin,Rout,Rhy;
  std::vector<std::vector<std::vector<point> > > allSegments;
  std::vector<std::vector<std::vector<int> > > allTruthTables;
  int solventsPerLink;
  double helixRatio;
  int ntrivs;
  std::vector<std::vector<point> > helixpts;
  std::vector<std::vector<point> > helPtsFlat;
  std::vector<std::vector<point> > solPtsFlat;
  std::vector<std::vector<double> > distanceSets;
  std::vector<std::pair<double,double> > distanceList;
  int solIntbound,molIntbound,solmolIntbound,nscat,nSize;
  double maxDistChange,kminVal,kmaxVal,maxDistVal,maxdist,maxScatDist,mutualDistCutOff;
  bool isTooClose;
};

#endif
