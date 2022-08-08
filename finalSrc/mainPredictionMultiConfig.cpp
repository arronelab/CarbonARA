#include "ktlMoleculeRandom.h"
#include "hydrationShellRandom.h"
#include "experimentalData.h"
#include <string.h> 


/***********************************************

  argv[1] scattering data file
  argv[2] sequence file location
  argv[3] initial prediction coord location  file 
  argv[4] paired distances file (can be empty)
  argv[5] fixed sections file (again can be empty)
  argv[6] number of structures 
  argv[7] request to apply hydrophobic covering WITHIN monomers will be a list of sections on which to apply it. Will say none if not.
  argv[8] request to apply hydrophobic covering BETWEEN monomers will be a list of pairs to try to hydropobically pair. Will say none if not.
  argv[9] kmin
  argv[10] kmax
  argv[11] Max number of fitting steps
  argv[12] prediction file 
  argv[13] scattering output file
**********************************************/


bool checkTransition(double &chiSqVal,double &chiSqCurr,double &uniformProb,int index,int &maxSteps){
  /*double tempFrac;
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
  }*/
  if(chiSqVal<chiSqCurr){
    return true;
  }else{
    return false;
  }
}

void generateAllDistances(ktlMolecule &molin,hydrationShellMinimal &hydrationShellin,std::vector<double> &molDistsin,std::vector<double> &solDistsin,std::vector<double> &solMolDistsin,std::vector<double> &meanHydrophilicDists,int &molSizein,int &noSolin,double &maxDistin){
  // calculate calpha-calpha distances
  molDistsin = molin.getDistSet();
  molSizein = molin.getNoAminos();
  std::sort(molDistsin.begin(),molDistsin.end());
  // get solvent molecules
  std::vector<std::vector<point> > solptsin = hydrationShellin.returnFlatSolList();
  // calculate solvent-solvent distances
  solDistsin = molin.getExactDistSet(solptsin);
  std::sort(solDistsin.begin(),solDistsin.end());
  std::cout<<"minimum  solvent distance "<<solDistsin[0]<<"\n";
  maxDistin = solDistsin[solDistsin.size()-1];
  noSolin=0;
  double solMax=4.0;
  double totalHydroOverlap=0.0;
  meanHydrophilicDists = molin.getHydrophobicDistance(solptsin,solMax);
  for(int i=0;i<meanHydrophilicDists.size();i++){
    totalHydroOverlap =totalHydroOverlap + meanHydrophilicDists[i];
  }
  std::cout<<"how much hydro "<<totalHydroOverlap<<"\n";
  for(int i=0;i<solptsin.size();i++){
    noSolin = noSolin + solptsin[i].size();
  }
  // calculate solvent-calpha distances
  solMolDistsin = molin.solMolDists(solptsin);
  std::sort(solMolDistsin.begin(),solMolDistsin.end());
}

double getOverlapPenalty(double &closestApproachDist,std::vector<double> &overlapDists){
  double distSumCurr=0.0;
  for(int l=0;l<overlapDists.size();l++){
    distSumCurr = distSumCurr + std::exp(std::abs(closestApproachDist-overlapDists[l]))-1.0;
  }
  if(overlapDists.size()>0){
    distSumCurr =0.02*(1.0/overlapDists.size())*distSumCurr;
  }
  return distSumCurr;
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
    // read in the sequence and seocndary structure
    std::stringstream ss;
    int ind1=i+1;
    ss<<ind1;
    const char* str = ss.str().c_str(); 
    ktlMolecule molTmp;
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
    strcat(coordinateLoc,str);
    strcat(coordinateLoc,".dat");
    molTmp.readInCoordinates(coordinateLoc);
    // fine the hydrophobic residues
    molTmp.getHydrophobicResidues();
    mol.push_back(molTmp);
  }

  
  bool doAll = false;
  /*********************************
   
   check which sections are being altered
   
   ********************************************/
  std::vector<std::vector<int> > fixedSecLists;
  std::ifstream fixedSecFile;
  // stuff inbetween
  for(int i=0;i<noStructures;i++){
    std::stringstream ss;
    int ind1=i+1;
    ss<<ind1;
    const char* str = ss.str().c_str(); 
    std::ifstream fixedSecFile;
    std::vector<int> fixedSecList;
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

     read in the scattering and set up the scattering model

   ******************************************/

  experimentalData ed(argv[1]);


  /*************************************************

    generate the hydration layer model
   
   *************************************************/
  // the hydration shell parameters
  double Rin= 6.0;
  double Rout=7.0;
  double RShell = 5.5;
  int ntrivs=6;
  double helixRatio=1.0;
  int solventsPerLink =1;

  // call the hydration shell code
  std::vector<hydrationShellMinimal> hydrationShells;
  for(int i=0;i<mol.size();i++){
    hydrationShellMinimal hydrationShell(mol[i],Rin,Rout,RShell,ntrivs,helixRatio,solventsPerLink,closestApproachDist,rmin,rmax,lmin);
     // generate the initial hydration shell
    hydrationShell.tubeParamList();
    hydrationShell.constructInitialState();
    hydrationShell.allOverlap();
    hydrationShells.push_back(hydrationShell);
    // load contact predictions or similar distance constraints
    mol[i].loadContactPredictions(argv[4]);
  }
 

  /*******************************************
    establish the initial fit
   ******************************************/

  // Generate the initial distance distribution, grab all distances
  
  // check internal overlap of sections and fill overlap list

  //list of all the distances
  std::vector<std::vector<double> > molDists;
  std::vector<std::vector<double> > solDists;
  std::vector<std::vector<double> > solMolDists;
  //
  std::vector<int> molSize;
  std::vector<int> noSol;
  
  // note we take the binsize based on the largest molecule

  // check for any overlaps
  std::vector<std::vector<double> > overlapDistSet;;
  std::vector<double> maxDistList; 
  for(int i=0;i<mol.size();i++){
    std::vector<double> overlapDists= mol[i].checkOverlapWithRad(closestApproachDist);
    overlapDistSet.push_back(overlapDists);
    // calculate all inter-molcular distances for the scattering model
    double maxDistTemp;
    std::vector<double>  molDistsTmp;
    std::vector<double> solDistsTmp;
    std::vector<double>  solMolDistsTmp;
    std::vector<double> meanHydrophilicDistsTmp;
    int molSizeTmp,noSolTmp;
    generateAllDistances(mol[i],hydrationShells[i],molDistsTmp,solDistsTmp,solMolDistsTmp,meanHydrophilicDistsTmp,molSizeTmp,noSolTmp,maxDistTemp);
    molDists.push_back(molDistsTmp);
    solDists.push_back(solDistsTmp);
    solMolDists.push_back(solMolDistsTmp);
    molSize.push_back(molSizeTmp);
    noSol.push_back(noSolTmp);
    maxDistList.push_back(maxDistTemp);
  }
   
  
  double maxDist = *std::max_element(std::begin(maxDistList), std::end(maxDistList));
  double  kmin = std::atof(argv[9]);
  double  kmax = std::atof(argv[10]);

  // set the number of bins using the shannon sampling theorem
  int noDistBins = int(1.1*std::ceil((kmax-kmin)*maxDist/3.14159265359));
  ed.setPhases(maxDist,kmin,kmax);

  // get the initial scattering

  std::vector<double> percentageCombinations(mol.size(),0.0);
  
  double scatterFit=10000.0;
  for(int i=8;i<=12;i++){
    std::vector<double> percentageCombinationsTmp(mol.size(),0.0);
    double p1 = double(i/20.0);
    double p2 = 1.0-p1;
    percentageCombinationsTmp[0]=p1;
    percentageCombinationsTmp[1]=p2;
    double scatterFitTemp = ed.fitToScatteringMultiple(molDists,solDists,solMolDists,molSize,noSol,percentageCombinationsTmp);
    std::cout<<p1<<" "<<p2<<" "<<scatterFitTemp<<"\n";
    if(scatterFitTemp<scatterFit){
      scatterFit = scatterFitTemp;
      percentageCombinations = percentageCombinationsTmp;
    }
  }
  //check for overlaps
  

  // Penalise excessive overlap

  double overlapPenalty = 0.0;

  for(int i=0;i<mol.size();i++){
    overlapPenalty = overlapPenalty + getOverlapPenalty(closestApproachDist,overlapDistSet[i]);
  }

  scatterFit =scatterFit + overlapPenalty;

  // add any additional constrainst to the overall fitness (should rename as its not necssarily Lennard-Jones

  double contactPredPen=0.0;
  for(int i=0;i<mol.size();i++){
    contactPredPen =contactPredPen + mol[i].getLennardJonesContact();
  }

  scatterFit = scatterFit +  contactPredPen;

  /******************************************************************
  
                      Fit to scattering data
  
  ******************************************************************/

  //ed.writeScatteringToFile(molDists,solDists,solMolDists,molSize,noSol,argv[13]);

  /****************************************************************************
    
    Main algorithm 
  
  ***************************************************************************/

  //set up loop parameters
  int k=0;
  std::vector<int> noSections;
  for(int i=0;i<mol.size();i++){
    int noSectionsTmp = mol[i].noChains();
    noSections.push_back(noSectionsTmp);
  }
  int noScatterFitSteps=std::atoi(argv[11]);
  //std::cout<<"no chains "<<noSections<<" "<<scatterFit<<" "<<noScatterFitSteps<<"\n";
  
  while(k<noScatterFitSteps && scatterFit>0.0001){
    // start random genrator (only necessary for simulated annealing)
    std::random_device rdev{};
    std::default_random_engine generator1{rdev()};
    std::uniform_real_distribution<double> distributionR(0.0,1.0);
    double scatterFitBest = scatterFit;
    // loop over the number of units of the molecule (i.e. 2 for a dimer one for a monomer)
    for(int l=0;l<mol.size();l++){
      int netIndex=0;
      for(int i=1;i<=noSections[l];i++){
	// net index tells us how far we are through the whole moelcule
	if(i>1){
	  netIndex=netIndex+mol[l].getSubsecSize(i-1);
	}
	// Now loop over the secondary structures of the given unit or section
	ktlMolecule molBest = mol[l];
	std::vector<std::vector<double> > molDistsBest = molDists;
	std::vector<std::vector<double> > solDistsBest = solDists;
	std::vector<std::vector<double> > solMolDistsBest = solMolDists;
	std::vector<int> noSolBest = noSol;
	double overlapPenaltyBest;
	int noDistBinsBest;
	std::vector<double> maxDistListBest = maxDistList;
	std::vector<double> percentageCombinationsBest(mol.size(),0.0);
	for(int j=0;j<mol[l].getSubsecSize(i)-1;j++){
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
	    // this checls if there haven't been any rouge sections created (some occasional flaws in the procedure which are to be ironed ou
	    bool cacaDist= molCopy.checkCalphas(i);
	    std::vector<std::vector<double> > molDistsTmp = molDists;
	    std::vector<std::vector<double> > solDistsTmp = solDists;
	    std::vector<std::vector<double> > solMolDistsTmp = solMolDists;
	    std::vector<int> noSolTmp = noSol;
	    if(cacaDist==false){
	      // if the new structure is okay we calculate its fit
	      // generate a new hydration shell
	      hydrationShellMinimal hydrationShellCopy(molCopy,Rin,Rout,RShell,ntrivs,helixRatio,solventsPerLink,closestApproachDist,rmin,rmax,lmin);
	      // generate the hydration shell
	      hydrationShellCopy.tubeParamList();
	      hydrationShellCopy.constructInitialState();
	      hydrationShellCopy.allOverlap();
	      std::vector<double> molDistsCopy;
	      std::vector<double> solDistsCopy;
	      std::vector<double> solMolDistsCopy;
	      std::vector<double> meanHydrophilicDistsTmp;
	      //
	      int molSizeTmp,noSolCopy;
	      double maxDistCopy;
	      // check for any overlaps
	      std::vector<double> overlapDistsNew= molCopy.checkOverlapWithRad(closestApproachDist);
	      generateAllDistances(molCopy,hydrationShellCopy,molDistsCopy,solDistsCopy,solMolDistsCopy,meanHydrophilicDistsTmp,molSizeTmp,noSolCopy,maxDistCopy);
	      std::vector<double> maxDistListTmp = maxDistList;
	      // very rarely the solvent produciton algo misbehaves so.  
	      if(solDistsCopy[solDistsCopy.size()-1]<10000.0){
		//new scattering
		int newNoDistBins;
		maxDistListTmp[l]=maxDistCopy;
		double maxDistCopy = *std::max_element(std::begin(maxDistListTmp), std::end(maxDistListTmp));
		newNoDistBins = int(1.1*std::ceil((kmax-kmin)*maxDistCopy/3.14159265359));
		//std::cout<<"bins "<<newNoDistBins<<" "<<noDistBins<<"\n";
		if(newNoDistBins!= noDistBins){
		  ed.setPhases(maxDist,kmin,kmax);
		}
		molDistsTmp[l] = molDistsCopy;
		solDistsTmp[l] = solDistsCopy;
		solMolDistsTmp[l] = solMolDistsCopy;
		noSolTmp[l] = noSolCopy;
		double newScatterFit=10000.0;
		 std::vector<double> percentageCombinationsTmp(mol.size(),0.0);
		for(int i=8;i<12;i++){
		   std::vector<double> percentageCombinationsTmpTmp(mol.size(),0.0);
		  double p1 = double(i/20.0);
		  double p2 = 1.0-p1;
		  percentageCombinationsTmpTmp[0]=p1;
		  percentageCombinationsTmpTmp[1]=p2;
		  double scatterFitTemp = ed.fitToScatteringMultiple(molDistsTmp,solDistsTmp,solMolDistsTmp,molSize,noSolTmp,percentageCombinationsTmpTmp);
		  //std::cout<<p1<<" "<<p2<<" "<<scatterFitTemp<<"\n";
		  if(scatterFitTemp<newScatterFit){
		    newScatterFit = scatterFitTemp;
		    percentageCombinationsTmp = percentageCombinationsTmpTmp;
		  }
		}
		double overlapPenaltyNew=0.0;
		double contactPredPenNew=0.0;
		for(int m=0;m<mol.size();m++){
		  if(m==l){
		    overlapPenaltyNew = overlapPenaltyNew + getOverlapPenalty(closestApproachDist,overlapDistsNew);
		    contactPredPenNew = contactPredPenNew + molCopy.getLennardJonesContact();
		  }else{
		    overlapPenaltyNew = overlapPenaltyNew + getOverlapPenalty(closestApproachDist,overlapDistSet[m]);
		    contactPredPenNew = contactPredPenNew + mol[m].getLennardJonesContact();
		  }
		}
		//std::cout<<" whats gwanin ? "<<l<<" "<<newScatterFit<<" "<<overlapPenaltyNew<<"\n";
		newScatterFit =newScatterFit +  overlapPenaltyNew;
		// finally add contract prediction value
		newScatterFit =newScatterFit + contactPredPenNew;
		
		double uProb = distributionR(generator1);
		if(checkTransition(newScatterFit,scatterFitBest,uProb,k,noScatterFitSteps)){
		  scatterFitBest  = newScatterFit;       
		  molBest=molCopy;
		  molDistsBest = molDistsTmp;
		  solDistsBest = solDistsTmp;
		  solMolDistsBest = solMolDistsTmp;
		  noSolBest = noSolTmp;
		  overlapPenaltyBest=overlapPenaltyNew;
		  noDistBinsBest=newNoDistBins;
		  maxDistListBest = maxDistListTmp;
		  percentageCombinationsBest = percentageCombinationsTmp;
		}
	      }
	    }
	  }
	}
	double uProb = distributionR(generator1);
	if(checkTransition(scatterFitBest,scatterFit,uProb,k,noScatterFitSteps)){
		scatterFit  = scatterFitBest;       
		mol[l]=molBest;
		molDists = molDistsBest;
		solDists = solDistsBest;
		solMolDists = solMolDistsBest;
		noSol = noSolBest;
		overlapPenalty=overlapPenaltyBest;
		noDistBins=noDistBinsBest;
		maxDistList = maxDistListBest;
		percentageCombinations = percentageCombinationsBest;
	}
      }
    }
    k++;
    std::cout<<k<<" "<<scatterFit<<" "<<overlapPenalty<<" "<<percentageCombinations[0]<<" "<<percentageCombinations[1]<<"\n";
  }
  for(int i=0;i<mol.size();i++){
    std::stringstream ss;
    int ind1=i+1;
    ss<<ind1;
    const char* str = ss.str().c_str(); 
    char outputMolLoc[100];
    strcpy(outputMolLoc,argv[12]);
    strcat(outputMolLoc,"_");
    strcat(outputMolLoc,str);
    strcat(outputMolLoc,".dat");
    mol[i].writeMoleculeToFile(outputMolLoc);
  }
  ed.writeScatteringToFileMultiple(molDists,solDists,solMolDists,molSize,noSol,percentageCombinations,argv[13]);
}
      
