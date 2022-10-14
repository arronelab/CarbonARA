/**
 * @file ktlMoleculeRandom.cpp
 *
 * @brief Primary Molecule Class: Reading, Writing, Manipulation, Generation
 *
 * @ingroup Generation
 *
 *
 * @author Chris Prior
 * Contact: christopher.prior@durham.ac.uk
 *
 */

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
 * Big con
 *
 * Tracing is controlled on a per "component" basis, where a "component" is a
 * name of the form aaa.bbb.ccc where aaa is the Most significant part; for
 * example, the utilities library might be called "utils", the doubly-linked
 * list "utils.dlist", and the code to destroy a list "utils.dlist.del"
 *
 */
class ktlMolecule{
public:
  ktlMolecule();
  //ktlMolecule(ktlMolecule &ktlm);
  void setParams(double &rminIn,double &rmaxIn,double &lminIn);
  void readInMol(const char* filename,int chainNo,int isRand);
  void readInMolGen(const char* filename);
  void readInMolWithBackbone(const char* filename,int chainNo,int isRand);
  void readInMolWithBackboneLenJ(const char* filename,int chainNo,int isRand,const char* LenJFileName);
  void readInMolWithSequenceLenJ(const char* filename,int chainNo,int isRand,const char* LenJFileName);
  std::vector<int> getUnitNos();
  std::vector<std::vector<point> > getTangents();
  std::vector<std::vector<point> > getNormals();
  std::vector<std::vector<point> > getBinormals();
  std::vector<std::vector<point> > getCoordinates();
  std::vector<point> getCoordinatesSection(int i);
/**
* Returns Size of Secondary Structure Element
*
*@param sec Index of Secondary Structure Element (Counting from 1)
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
* Returns Number of Chains
*/
  int noChains();
  int noSecSize();
  int getNoAminos();
  void createSyntheticMolecule(int nLinks);
  void readInFit(const char* filename,const char* fieldname);
  void readInMolNmer(const char* filename,std::vector<int> &molIndicies,int isRand,double &rmin,double &rmax,double &lmin);
  void readInMolNmerAndCoordinates(const char* filename,std::vector<int> &molIndicies,int isRand,double &rmin,double &rmax,double &lmin);
  void readInMolNmerWithBackbone(const char* filename,std::vector<int> &molIndicies,int isRand,double &rmin,double &rmax,double &lmin);
  void readInMolNmerWithBackboneLenJ(const char* filename,std::vector<int> &molIndicies,int isRand,double &rmin,double &rmax,double &lmin,const char* filenameLenJ);
  void readInMolNmerWithSequenceLenJ(const char* filename,std::vector<int> &molIndicies,int isRand,double &rmin,double &rmax,double &lmin,const char* filenameLenJ);
  /** 
* Processes Sequence File
* 
* Reads in Number of Chains, Primary Sequence, and (Simplified) Secondary Sequence. Uses Secondary Sequence to partition the curve.
*
*@param filename Location of Sequence Fingerprint file 
*
*@param rmin Minimum distance between neighbouring cAlpha residues, default=3.7
*
*@param rmax  Maximum distance between neighbouring cAlpha residues, default=3.9
*
*@param lmin Closest distance two non adjactent local (same secondary structure) cAlpha residues can get, default=4.0
*/
  void readInSequence(const char* filename,double &rmin,double &rmax,double &lmin);
  /**
* Processes Coordinate File
*
* Description
*
*@param filename Location of Coordinates file
*/
  void readInCoordinates(const char* filename);
  void readInSequenceWBackbone(const char* filename,int chainNo,const char* backbonename);
  // void readInMolNmerSequence(const char* filename,std::vector<int> &molIndicies,double &rmin,double &rmax,double &lmin);
  void readInMolNmerSequenceWBackbone(const char* filename,const char* backboneName,std::vector<int> &molIndicies,double &rmin,double &rmax,double &lmin);
  point binorm(point &T1,point &T2);
  point parallelTransport(point &tan1,point &tan2,point &norm1);
  std::vector<std::vector<point> > updateFrame(std::vector<point> &section,point &tangent,point &normal,point &binormal);
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
  double getGlobalRadiusOfCurvature();
  double getGlobalRadiusOfCurvatureWithinSec(int &secNo,int &NoNeighbour);
  double getGlobalRadiusOfCurvatureBetweenSec(int &NoNeighbour);
  point getCentreOfMass(std::vector<std::vector<point> > &cdSet);
  std::vector<int> checkOverlap(std::vector<std::vector<point> > &cdsIN);
  std::vector<double> checkOverlapWithRad(double &wRad,int &sec);
  std::vector<double> checkOverlapWithRad(double &wRad);
  std::vector<double> getDistSet();
  double compareDistances(std::vector<std::vector<point> > &coords2);
/**
*
*
*
*/
  bool checkCalphas(std::vector<std::vector<point> > &coordsIn);
  bool checkCalphas();
  bool checkCalphas(int &index);
  int getRandomMolecule();
  void getRandomMoleculeAllowOverlap();
  int getRandomMoleculeReset();
  void getRandomMoleculeAllowOverlapReset();
  void removeOverlap();
  void resetRandomMolecule();
  void getFrameForBackbone();

/** 
* Resampling of single molecule section
* 
* Some description of how that is done...
*
* @param index Index of protein section to be re-generated
*/
  void changeMoleculeSingle(int &index);


  void changeMoleculeSet(std::vector<int> &indicies);
/** 
* Changes singe secondary structure element
* 
* Some description of how that is done...
*
*@param index Index of secondary structure element to be changed
*
*@param sec Chain number
*/
  void changeMoleculeSingleMulti(int &index,int sec);
  void changeMoleculeSetMulti(std::vector<int>  &indicies,int sec);
  void changeMoleculeMultiRotate(double &angle,point &k,int secIn,point &transVec);
  void replicateMolecule(int &noReplications); 
  void changeMoleculeLocal(int &index,double variationSize);
  void changeMoleculeLocalSet(std::vector<int> &indicies,double variationSize);
  void changeMoleculeSingle(int &index,std::vector<std::vector<point> > &cdsIn,std::vector<std::pair<std::string,int> > &nameSizeSubList);
  void rotation3D(point &p,point &centre,point &k,double &cosangle,double &sinangle);
  void rotateSection(std::vector<std::vector<point> >  &section,point &centre,point &k,double &angle,point &transVec);
  bool checkOverlapSbond(double &minDist);
  void changeMoleculeSingleCheckOverlap();
/**
* Outputs molecule as a coordinate file
*
* Multimers will be in one single file, with each chain separate by "End Chain i"
*
* @param filename Ouput coordinate file location
*/


  void writeMoleculeToFile(const char* filename);
  void getBackboneStats();
  std::vector<std::pair<double,double> > getKapTauVals();
  double getPairDistance(std::pair<int,int> &index1,std::pair<int,int> &index2);
  double lennardJones(double &idealDist,double &currDist,int noConPred,double &weightCoeff);
  double getFitQualSbonds(std::vector<std::tuple<std::pair<int,int>,std::pair<int,int>,double> > &contactPairList,double weightCoeff);
  double getFitQualSbondsNonSpecific(std::vector<std::pair<std::pair<int,int>,double >> &contactPairList,double weightCoeff);
  std::vector<double> getExactDistSet(std::vector<std::vector<point> > &coordSet);
  double allDists(std::vector<point> &ptSet,std::vector<double> &distSet);
  std::vector<double> solMolDists(std::vector<std::vector<point> > &pts1);
  std::vector<double> getApproxDistSet(double &divFac);
  std::vector<double> getApproxDistSetFixed();
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
