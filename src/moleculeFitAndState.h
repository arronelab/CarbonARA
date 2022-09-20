#ifndef ALG_ROUTE
#define ALG_ROUTE

#include "ktlMoleculeRandom.h"
#include "hydrationShellRandom.h"
#include "localWrithe.h"
#include "experimentalData.h"

class moleculeFitAndState{
public:
  moleculeFitAndState(std::vector<ktlMolecule> &mol,double &RinSv,double &RoutSv,double &RShellSv,int &ntrivsSv,double &closestApproachDistIn,int &solventsPerLinkIn,double &rminIn,double &rmaxIn,double &lminIn);
  std::vector<ktlMolecule> getMolecule();
  void updateMolecule(std::vector<ktlMolecule> &molNew);
  double calculateScattering(experimentalData &ed,double &kmin,double &kmax,std::vector<double> &mixtureVals);
  void writeScatteringToFile(experimentalData &ed,double &kmin,double &kmax,const char* filename);
  double getOverlapPenalty(double &closestApproachDist,std::vector<double> &overlapDists);
  double applyOverlapPenalty();
  double applyDistanceConstraints();
  double applyDistanceConstraints(ktlMolecule &molNew,int &i);
  void calculateMoleculeDistances(ktlMolecule &molNew,int &i);
  void calcuateHydrationDistances(hydrationShellMinimal &hs,int &i);
  void applyWritheConstraint();
  double getFit();
  void alterWritheSet(ktlMolecule &molNew,int &i);
  double getOverallFit(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,std::vector<double> &helRatList,double &kmin,double &kmax);
  double getOverallFit(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,std::vector<double> &helRatList,ktlMolecule &molNew,double &kmin,double &kmax,int &i);
  double currFit;
private:
  std::vector<std::vector<double> > molDists;
  std::vector<std::vector<double> > solDists;
  std::vector<std::vector<double> > solMolDists;
  std::vector<std::vector<double> > overlapDistSet;
  std::vector<int> molSize;
  std::vector<int> noSol;
  double maxDist;
  double hydroPhobicPacking;
  std::vector<double> originalWrithes;
  std::vector<double> currWrithes;
  std::vector<double> maxDistMol;
  std::vector<double> maxDistSol;
  std::vector<double> contactPredPen;
  double writhePenalty;
  double Rin,Rout,RShell,ntrivs,closestApproachDist;
  double solventsPerLink,rmin,rmax,lmin;
  std::vector<double> percentageCombinations;
  std::vector<ktlMolecule> mol;
};

#endif
