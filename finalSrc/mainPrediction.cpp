#include "ktlMoleculeRandom.h"
#include "hydrationShellRandom.h"
#include "experimentalData.h"
#include <string.h>
#include "localWrithe.h"
#include "moleculeFitAndState.h"


/***********************************************

  argv[1] scattering data file
  argv[2] sequence file location
  argv[3] initial prediction coord location  file 
  argv[4] paired distances file (can be empty)
  argv[5] fixed sections file (again can be empty)
  argv[6] number of structures 
  argv[7] request to apply hydrophobic covering WITHIN monomers will be a list of sections on which to apply it. Will say none if not. -- Currently not used
  argv[8] request to apply hydrophobic covering BETWEEN monomers will be a list of pairs to try to hydropobically pair. Will say none if not. -- currently not used
  argv[9] kmin 
  argv[10] kmax
  argv[11] Max number of fitting steps
  argv[12] prediction file 
  argv[13] scattering output file
  argv[14] mixture list file, a list of sets of numbers indicatig the allowed set of mixture percentages of each species (e.g. dimer 20 monomer 80)
  argv[15] target water ratio
**********************************************/


bool checkTransition(double &chiSqVal,double &chiSqCurr,double &uniformProb,int index,int &maxSteps){
  double tempFrac;
  if(index < maxSteps/2){
    tempFrac= 1.0-double(index)/double(maxSteps/2);
  }else{
    tempFrac=0.0000000001;
  }
  double annealProb;
  double tempVal;
  if(chiSqVal<chiSqCurr){
    annealProb=1.0;
  }else{
    annealProb = std::exp(-(chiSqVal-chiSqCurr)/tempFrac);
  }
  //std::cout<<index<<" "<<maxSteps<<" "<<chiSqVal<<" "<<chiSqCurr<<" "<<annealProb<<" "<<uniformProb<<"\n";
  if(annealProb>uniformProb){
    return true;
  }else{
    return false;
  }
  /* if(chiSqVal<chiSqCurr){
    return true;
  }else{
    return false;
    }*/
}

void sortVec(std::vector<moleculeFitAndState> &mfs){
  std::sort(mfs.begin(), mfs.end(),[](const moleculeFitAndState &x, const moleculeFitAndState &y) {
                
    return x.currFit < y.currFit;
  });
}



double getHydrophobicPackingPenalty(double &packValue){
  return 0.00001*std::exp(3.0*(packValue-1.6));
}

int main( int argc, const char* argv[] )
{  
 /*************************************

  set up model parameters 
 
  *************************************/
  
  double lmin=4.0; // closest distance two non adjactent local (same secondary unit) moelcules can get
  double rmin=3.7;double rmax=3.9; // max and min calpha-Dists
  double closestApproachDist=3.9; // closest distance two non adjactent non local moelcules (different secondary unit) can get

  /*************************************
  
   determine initial model: Two options no initial prediction, we must generate a structure
   or some initial structure provided. Actually we need a half-half option

   *************************************/

  int noStructures = std::atoi(argv[6]);
  
  
  std::vector<ktlMolecule> mol;
  for(int i=0;i<noStructures;i++){
    // read in the sequence and secondary structure
    std::stringstream ss;
    ktlMolecule molTmp;
    int ind1=i+1;
    ss<<ind1;
    //std::cout<<ind1<<"\n";
    const char* str = ss.str().c_str();
    char sequenceLoc[100];
    strcpy(sequenceLoc,argv[2]);
    strcat(sequenceLoc,"fingerPrint");
    strcat(sequenceLoc,str);
    strcat(sequenceLoc,".dat");
    molTmp.readInSequence(sequenceLoc,rmin,rmax,lmin);
    // read in the initial coordinates
    char coordinateLoc[100];
    strcpy(coordinateLoc,argv[2]);
    strcat(coordinateLoc,"coordinates");
    str = ss.str().c_str();
    strcat(coordinateLoc,str);
    strcat(coordinateLoc,".dat");
    molTmp.readInCoordinates(coordinateLoc);
    // fine the hydrophobic residues
    molTmp.getHydrophobicResidues();
    mol.push_back(molTmp);
  }

  
  bool doAll = false;
  
  
  /*********************************
   
   Determine which sections are being altered
   
   ********************************************/
  std::vector<std::vector<int> > fixedSecLists;
  std::ifstream fixedSecFile;
  // stuff inbetween
  for(int i=0;i<noStructures;i++){
    std::stringstream ss;
    std::ifstream fixedSecFile;
    std::vector<int> fixedSecList;
    int ind1=i+1;
    ss<<ind1;
    const char* str = ss.str().c_str(); 
    char fixedSecLoc[100];
    strcpy(fixedSecLoc,argv[2]);
    strcat(fixedSecLoc,"varyingSectionSecondary");
    strcat(fixedSecLoc,str);
    strcat(fixedSecLoc,".dat");
    fixedSecFile.open(fixedSecLoc);
    std::string line;int index;
    if(fixedSecFile.is_open()){
      while(!fixedSecFile.eof()){
	std::getline(fixedSecFile,line);
	std::stringstream ss(line);
	ss>>index;
	fixedSecList.push_back(index);
      }
    }else{
      std::cout<<"failed to open fixed section file\n";
    }
    fixedSecLists.push_back(fixedSecList);
    fixedSecFile.close();
  }
  
  /******************************************

     Read in the scattering and set up the scattering model

   ******************************************/

  experimentalData ed(argv[1]);


  /******************************************

   Read in the permissible mixture list
 
   *********************************************/ 

  std::vector<std::vector<double> > mixtureList;

  // read the allowed mixtures in from file
  std::ifstream permissibleMixtureFiles;
  permissibleMixtureFiles.open(argv[14]);
  std::string line;double index;
  if(permissibleMixtureFiles.is_open()){
    while(!permissibleMixtureFiles.eof()){
      std::vector<double> mixtureSet;
      std::getline(permissibleMixtureFiles,line);
      std::stringstream ss(line);
      //std::cout<<line<<" "<<line.length()<<"\n";
      while(ss.eof()==false){
	ss>>index;
	mixtureSet.push_back(index);
      }
      if(line.length()>0){
	mixtureList.push_back(mixtureSet);
      }
    }
  }else{
      std::cout<<"failed to open mixture file\n";
  }
  
  

  
  /*************************************************

  Initialise  hydration shell parameters 
   
   *************************************************/
  // the hydration shell parameters
  double Rin= 6.0;
  double Rout=7.0;
  double RShell = 5.5;
  int ntrivs=6;
  int solventsPerLink =1;
  std::vector<double> helRatList;
  helRatList.push_back(0.5);
  // scattering min max values
  double  kmin = std::atof(argv[9]);
  double  kmax = std::atof(argv[10]);
  
  // calculate the initial fit, this includes all constraints

  /*************************************************

   initialise the first fit 
   
   *************************************************/

  moleculeFitAndState molFit(mol,Rin,Rout,RShell,ntrivs,closestApproachDist,solventsPerLink,rmin,rmax,lmin);
  double scatterFit= molFit.getOverallFit(ed,mixtureList,helRatList,kmin,kmax);
 
  /****************************************************************************
    
     Main algorithm 
  
  ***************************************************************************/
  // 


  // //set up loop parameters
  int k=0;


  
  // /* the vector noSections tells us how many subsections are in each moelcule
  //    e.g. for a monomer/dimer mixture noSections[0]=1,noSections[1]=2.
  //  */
  
  std::vector<int> noSections;
  for(int i=0;i<mol.size();i++){
    int noSectionsTmp = mol[i].noChains();
    noSections.push_back(noSectionsTmp);
  }
  
  int noScatterFitSteps=std::atoi(argv[11]);


  // set up the vector of existing structures

  std::vector<moleculeFitAndState> molFitAndStateSet;

  // initialise to the initial guess
  int noMoleculeHist =10;
  for(int i=0;i<noMoleculeHist;i++){
    molFitAndStateSet.push_back(molFit);
  }
  
  std::random_device rdev{};
  std::default_random_engine generator{rdev()};
  while(k<noScatterFitSteps && scatterFit>0.00001){

    // loop over the molecules (e.g monomer and dimer fit

    std::uniform_real_distribution<double> distributionR(0.0,1.0);
    
    //choose which fit to
    double p = 0.7-0.6*(k/noScatterFitSteps);
    std::binomial_distribution<> changeIndexProbability(noMoleculeHist-1,p);
    int index =  changeIndexProbability(generator);
    molFit = molFitAndStateSet[index];
    //std::cout<<"checking mol no "<<index<<"\n";
    mol = molFit.getMolecule();
    scatterFit = molFit.currFit;
    for(int l=0;l<mol.size();l++){
      int netIndex=0;

      //loop over the sections of the given molecule (i.e. if its a monomer this loop is tivial, but not for a multimer
      
      for(int i=1;i<=noSections[l];i++){
	
	// net index tells us how far we are through the whole moelcule
	if(i>1){
	  netIndex=netIndex+mol[l].getSubsecSize(i-1);
	}
	
	// Now loop over the secondary structures of the given unit or section

	for(int j=0;j<mol[l].getSubsecSize(i)-1;j++){
	  //std::cout<<" mol "<<l<<" sec "<<i<<" has this many sections "<<mol[l].getSubsecSize(i)<<"\n";
	  int totalIndex = netIndex+j;
	  // in this if statement we check which secondary sections are being changed 
	  if((doAll==true) || (std::find(fixedSecLists[l].begin(),fixedSecLists[l].end(),totalIndex)!=fixedSecLists[l].end())){
	    // print statement currently in to check what we are changing is correct
	    //std::cout<<" section "<<totalIndex<<" of unit "<<i<<" "<<" sub set number "<<totalIndex-netIndex<<" being altered "<<mol.getSubsecSize(i)<<"\n";
	    
	    
	    // copy the molecule to change it and test if we do better
	    
	    ktlMolecule molCopy = mol[l];
	    
	    // now we change the section
	    int indexCh = totalIndex-netIndex;
	    molCopy.changeMoleculeSingleMulti(indexCh,i);
	    
	    // this (checkCalphas) checks if there haven't been any rouge sections created (some occasional flaws in the procedure which are to be ironed out
	    bool cacaDist= molCopy.checkCalphas(i);
	    if(cacaDist==false){

	      // calculate the new fit for this
	      moleculeFitAndState molFitTmp = molFit ;
	      //calculate all amino acid distances for changed molecule
	      double fitTemp = molFitTmp.getOverallFit(ed,mixtureList,helRatList,molCopy,kmin,kmax,l);
	      // check if we have imporved
	      // std::cout<<fitTemp<<" "<<scatterFit<<"\n";
	      double uProb = distributionR(generator);
	      if(checkTransition(fitTemp,scatterFit,uProb,k,noScatterFitSteps)){
		scatterFit = fitTemp;
		mol[l]=molCopy;
		molFit= molFitTmp;
	      }
	    }
	  }
	}
      }
    }
    molFitAndStateSet[index] = molFit;
    molFitAndStateSet[index].updateMolecule(mol);
    sortVec(molFitAndStateSet);
    for(int i=0;i<noMoleculeHist;i++){
      std::cout<<"step "<<k<<" "<<i<<" "<<molFitAndStateSet[i].currFit<<"\n";
    }
    k++;
  }
  mol =molFitAndStateSet[0].getMolecule();
  for(int i=0;i<mol.size();i++){
    std::stringstream ss;
    int ind1=i+1;
    ss<<ind1;
    char outputMolLoc[100];
    strcpy(outputMolLoc,argv[12]);
    strcat(outputMolLoc,"_");
    const char* str = ss.str().c_str(); 
    strcat(outputMolLoc,str);
    strcat(outputMolLoc,".dat");
    mol[i].writeMoleculeToFile(outputMolLoc);
  }
  // regenrate molecule hydration layer to update tha fit
  moleculeFitAndState molFitOut(mol,Rin,Rout,RShell,ntrivs,closestApproachDist,solventsPerLink,rmin,rmax,lmin);
  double scatterFitOut= molFitOut.getOverallFit(ed,mixtureList,helRatList,kmin,kmax);
  molFitOut.writeScatteringToFile(ed,kmin,kmax,argv[13]);
}
      
