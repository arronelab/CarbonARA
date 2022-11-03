#ifndef KTL_MOL
#define KTL_MOL

#include "point.h"
#include <fstream>
#include <sstream>
#include <string.h>
#include <algorithm>
#include <dirent.h>
#include "polyHelix.h"
#include <random>
#include "randomMolGen.h"
#include <tuple>

 /**
 * 
 * Kappa - Tau - Length - Molecule
 *
 * Do you have a molecule? \n
 * Well this is the class for you! \n
 * Read in the a molecule's coordinates and meta data. \n
 * Find geometric properties of your molecule's CA backbone. \n
 * Manipulate & DOMINATE your molecule, bend it to your will. \n
 * It will fit that SAXS curve... one way or another. \n
 * 
 */
class ktlMolecule{
public:
  ktlMolecule();
  //ktlMolecule(ktlMolecule &ktlm);
  void setParams(double &rminIn,double &rmaxIn,double &lminIn);
  std::vector<int> getUnitNos();
  std::vector<std::vector<point> > getTangents();
  std::vector<std::vector<point> > getNormals();
  std::vector<std::vector<point> > getBinormals();
  std::vector<std::vector<point> > getCoordinates();
  std::vector<point> getCoordinatesSection(int i);

  /**
  * Returns Size of Secondary Structure Element
  *
  * @param sec Index of Secondary Structure Element (Counting from 1)
  */ 
  int getSubsecSize(int sec);
  std::vector<std::vector<point> > getSubsecCoordinates(int &sec);
  std::vector<std::pair<std::string,int> > getNameSizeListOfSection(int &sec);
  //std::vector<point> getEulerAngles();
  std::vector<double> getDistChanges();
  double getCurvatureJoined(int index);
  double getTorsionJoined(int index);
  double getLengthJoined(int index);
  int getUnitNo(int index);
  double maxNeighbourDistSec(int &sec);
  point getTangent(int mindex,int subindex);
  point getNormal(int mindex,int subindex);
  point getBinormal(int mindex,int subindex);
  point getCoordinate(int mindex,int subindex);
  double getAlbeadJoined(int index);
  std::string getType(int &index);
  std::string getType(int &chainNo,int &index);
  //double getDistChange(int index);
  double getMaxDistChange();
  std::pair<double,double> getMaxPossibleLength();
  int molSize();

  /**
  * Returns Number of Chains (unexpected?)
  */
  int noChains();
  int noSecSize();
  int getNoAminos();
  void createSyntheticMolecule(int nLinks);

  /** 
  * Processes Sequence File
  * 
  * Reads in Number of Chains, Primary Sequence, and (Simplified) Secondary Sequence. Uses Secondary Sequence to partition the curve.
  *
  * @param filename Location of Sequence Fingerprint file 
  * @param rmin Minimum distance between neighbouring cAlpha residues, default=3.7
  * @param rmax  Maximum distance between neighbouring cAlpha residues, default=3.9
  * @param lmin Closest distance two non adjactent local (same secondary structure) cAlpha residues can get, default=4.0
  * 
  */
  void readInSequence(const char* filename,double &rmin,double &rmax,double &lmin);

  /**
  * Processes Coordinate File
  *
  *@param filename Location of Coordinates file
  */
  void readInCoordinates(const char* filename);

  /**
  * Finds Hydrophobic Residues
  *
  * Stores the index of all Hydrophobic Residues (Amino Acids: A,I,L,M,F,V,P,G)
  */
  void getHydrophobicResidues();
  std::vector<double> getHydrophobicDistance(std::vector<std::vector<point> > &solventList,double &maxSolDist);
  void getCoiledCoilResidues();
  void getPositiveResidues();
  void getNegativeResidues();
  void getPolarResidues();
  double coiledCoilPotential();
  double coiledCoilPotentialBetween(int &secNo);
  double coiledCoilPotentialBetween();
  point getCentreOfMass(std::vector<std::vector<point> > &cdSet);

  std::vector<int> checkOverlap(std::vector<std::vector<point> > &cdsIN);
  std::vector<double> checkOverlapWithRad(double &wRad);
  std::vector<double> checkOverlapWithRad(double &wRad,int &sec);

  std::vector<double> getDistSet();
  double compareDistances(std::vector<std::vector<point> > &coords2);

  bool checkCalphas();
  bool checkCalphas(int &index);
  bool checkCalphas(std::vector<std::vector<point> > &coordsIn);

  /** 
  * Resampling of single molecule section
  * 
  * Some description of how that is done...
  *
  * @param index Index of protein section to be re-generated
  * @param ?matrix of points?
  * @param need to check whats changed here
  * 
  */
  void changeMoleculeSingle(int &index,std::vector<std::vector<point> > &cdsIn,std::vector<std::pair<std::string,int> > &nameSizeSubList);

  int getRandomMolecule();
  int getRandomMoleculeReset();
  void resetRandomMolecule();

  /** 
  * Changes single secondary structure element
  * 
  * Some description of how that is done...
  *
  * @param index Index of secondary structure element to be changed
  * @param sec Chain number
  * 
  */
  void changeMoleculeSingleMulti(int &index,int sec);
  void changeMoleculeMultiRotate(double &angle,point &k,int secIn,point &transVec);
  void replicateMolecule(int &noReplications); 
  void rotation3D(point &p,point &centre,point &k,double &cosangle,double &sinangle);
  void rotateSection(std::vector<std::vector<point> >  &section,point &centre,point &k,double &angle,point &transVec);

  /**
  * Outputs molecule as a coordinate file
  *
  * Multimers will be in one single file, with each chain separate by "End Chain i"
  *
  * @param filename Ouput coordinate file location
  */
  void writeMoleculeToFile(const char* filename);
  std::vector<std::pair<double,double> > getKapTauVals();
  double getPairDistance(std::pair<int,int> &index1,std::pair<int,int> &index2);
  double lennardJones(double &idealDist,double &currDist,int noConPred,double &weightCoeff);
  std::vector<double> solMolDists(std::vector<std::vector<point> > &pts1);
  void loadContactPredictions(const char* contactloc);
  double getLennardJonesContact();
  void loadFixedSections(const char* fixedsecloc);
private:
  std::vector<int> noPts;
  std::vector<std::vector<point> > coords;
  std::vector<std::vector<std::string> > aminoList; 
  std::vector<std::vector<point> > tanlist;
  std::vector<std::vector<point> > normlist;
  std::vector<std::vector<point> > binormlist;
  polyHelix ph;
  std::vector<std::pair<int,int> > chainList;
  std::vector<double> distChanges;
  std::vector<std::pair<std::string,int> > nameSizeList;
  std::vector<std::pair<int,int> > hydroPhobicList;
  std::vector<std::pair<int,int> > polarList;
  std::vector<std::pair<int,int> > posChargeList;
  std::vector<std::pair<int,int> > negChargeList;
  std::vector<std::pair<int,int> > coiledCoilList;
  randomMol rmg;
  std::vector<std::vector<double> > distSetsSecs;
  std::vector<double> distSets;
  double kapvallink,kapvalbeta,kapvalalpha,tauvallink,tauvalbeta,tauvalalpha,alvallink,alvalbeta,alvalalpha,maxDistChange;
  std::vector<std::tuple<std::pair<int,int>,std::pair<int,int>,double> > contactPairList;
  std::vector<int> unchangedSections;
};

#endif
