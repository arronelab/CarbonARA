#ifndef ALG_ROUTE
#define ALG_ROUTE

#include "ktlMoleculeRandom.h"
#include "hydrationShellRandom.h"
#include "writhe.h"
#include "experimentalData.h"


/**
 * 
 * For checking structural metrics; proximity, writhe and fit to SAXS data
 * 
 */
class moleculeFitAndState{
public:
  moleculeFitAndState(std::vector<ktlMolecule> &mol,double &RinSv,double &RoutSv,double &RShellSv,int &ntrivsSv,double &closestApproachDistIn,int &solventsPerLinkIn,double &rminIn,double &rmaxIn,double &lminIn);
/**
* Calls on the mol class from ktlMolecule.
* 
* @return std::vector<ktlMolecule> 
*/
  std::vector<ktlMolecule> getMolecule();
/**
 * Updates the molecule to the current best fit, via ktlMolecule.
 * 
 * @param molNew The current best fit molecule.
 */
  void updateMolecule(std::vector<ktlMolecule> &molNew);
  double calculateScattering(experimentalData &ed,double &kmin,double &kmax,std::vector<double> &mixtureVals);
/**
 * Writes the scattering profile to given output file.
 * 
 * @param ed Inherits the experimentalData class, which handles scattering calculations.
 * @param kmin Minimum q value of the scattering data.
 * @param kmax Maximimum q value of the scattering data.
 * @param filename Location of output file.
 */
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
/**
 * Returns the goodness of fit for the current molecule, taking into account:
 * scattering data
 * a strict overlap penalty
 * an optional contact distance constraint penalty
 * a topological penalty which ensures the protein remains realistically entangled
 * 
 * @param ed Inherits the experimentalData class, which handles scattering calculations.
 * @param mixtureList Describes the relative ratios of monomer to multimer in the solution. If just looking at a monomer, will be simply 1.
 * @param helRatList No idea
 * @param kmin Minimum q value of the scattering data.
 * @param kmax Maximimum q value of the scattering data.
 * @return double 
 */
  double getOverallFit(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,std::vector<double> &helRatList,double &kmin,double &kmax);
  double getOverallFit(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,std::vector<double> &helRatList,ktlMolecule &molNew,double &kmin,double &kmax,int &i);
/**
 * Tracks the current fit to the scattering data
 * 
 */
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
