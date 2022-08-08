#include "ktlMoleculeRandom.h"
#include "hydrationShellRandom.h"
#include "experimentalData.h"
#include <string.h>
#include "localWrithe.h"


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

void generateAllDistances(ktlMolecule &molin,hydrationShellMinimal &hydrationShellin,std::vector<double> &molDistsin,std::vector<double> &solDistsin,std::vector<double> &solMolDistsin,std::vector<double> &meanHydrophilicDists,int &molSizein,int &noSolin,double &maxDistin,double &hydroPhobicPacking){
  // calculate calpha-calpha distances
  molDistsin = molin.getDistSet();
  molSizein = molin.getNoAminos();
  std::sort(molDistsin.begin(),molDistsin.end());
  // get solvent molecules
  std::vector<std::vector<point> > solptsin = hydrationShellin.returnFlatSolList();
  // calculate solvent-solvent distances
  std::sort(solDistsin.begin(),solDistsin.end());
  maxDistin = solDistsin[solDistsin.size()-1];
  noSolin=0;
  double solMax=5.5;
  hydroPhobicPacking =0.0;
  meanHydrophilicDists = molin.getHydrophobicDistance(solptsin,solMax);
  for(int i=0;i<meanHydrophilicDists.size();i++){
     hydroPhobicPacking  =  hydroPhobicPacking  + meanHydrophilicDists[i];
  }
  //std::cout<<"how much hydro "<<totalHydroOverlap<<"\n";
  for(int i=0;i<solptsin.size();i++){
    noSolin = noSolin + solptsin[i].size();
  }
  //std::cout<<"no amino acis "<<molSizein<<" no solvents "<<noSolin<<"\n";
  // calculate solvent-calpha distances
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
   
   Determine which sections are being altered
   
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

     Read in the scattering and set up the scattering model

   ******************************************/

  experimentalData ed(argv[1]);


  /*************************************************

    generate the hydration layer model(s)
   
   *************************************************/
  // the hydration shell parameters
  double Rin= 6.0;
  double Rout=7.0;
  double RShell = 5.5;
  int ntrivs=6;
  double helixRatio=0.5;
  int solventsPerLink =1;

  //list of all the distances
  std::vector<std::vector<double> > molDists;
  std::vector<std::vector<double> > solDists(mol.size());
  std::vector<std::vector<double> > solMolDists(mol.size());
  std::vector<int> molSize;
  std::vector<int> noSol;

  // call the hydration shell for eaxm molecule 
  std::vector<hydrationShellMinimal> hydrationShells;
  for(int i=0;i<mol.size();i++){
    hydrationShellMinimal hydrationShell(mol[i],Rin,Rout,RShell,ntrivs,helixRatio,solventsPerLink,closestApproachDist,rmin,rmax,lmin);
     // generate the initial hydration shell
    hydrationShell.tubeParamList();
    hydrationShell.constructInitialState();
    hydrationShell.getAllHelices();
    // Calculate the solvent to molecule distances 
    hydrationShell.solventMoleculeDistances(solMolDists[i],solDists[i]);
    hydrationShells.push_back(hydrationShell);
    // load contact predictions or similar distance constraints
    mol[i].loadContactPredictions(argv[4]);
  }
 

  /*******************************************
    establish the initial fit
   ******************************************/

  // Generate the initial distance distribution, grab all distances


  // overlapDistSet stores a list of distances below lmin, the smallest non-local distance
  std::vector<std::vector<double> > overlapDistSet;

  // note we take the binsize based on the largest molecule distance over all moelcules. maxDistList stores this
  std::vector<double> maxDistList;

  // The vector hydrophobicPackingMeasures stores hhyrdophoic packingscores for each moelcule (the exact nature of the score is work in progress)
  std::vector<double> hydrophobicPackingMeasures;

  // the main distance calculating loop,
  for(int i=0;i<mol.size();i++){
    //check for overlaps
    std::vector<double> overlapDists= mol[i].checkOverlapWithRad(closestApproachDist);
    overlapDistSet.push_back(overlapDists);
    // calculate all inter-molcular distances for the scattering model
    double maxDistTemp;
    std::vector<double>  molDistsTmp;
    //std::vector<double> solDistsTmp;
    //std::vector<double>  solMolDistsTmp;
    std::vector<double> meanHydrophilicDistsTmp;
    int molSizeTmp,noSolTmp;
    //std::cout<<i<<" "<<" here? \n";
    double hydroPhobicPackingMeasure;
    generateAllDistances(mol[i],hydrationShells[i],molDistsTmp,solDists[i],solMolDists[i],meanHydrophilicDistsTmp,molSizeTmp,noSolTmp,maxDistTemp,hydroPhobicPackingMeasure);
    molDists.push_back(molDistsTmp);
    molSize.push_back(molSizeTmp);
    noSol.push_back(noSolTmp);
    maxDistList.push_back(maxDistTemp);
    hydrophobicPackingMeasures.push_back(hydroPhobicPackingMeasure);
  }
   
  //fine the biggets distance to slect bing size
  double maxDist = *std::max_element(std::begin(maxDistList), std::end(maxDistList));
  double  kmin = std::atof(argv[9]);
  double  kmax = std::atof(argv[10]);

  // set the number of bins using the shannon sampling theorem
  int noDistBins = int(1.1*std::ceil((kmax-kmin)*maxDist/3.14159265359));
  ed.setPhases(maxDist,kmin,kmax);

  // get the initial scattering

  std::vector<double> percentageCombinations(mol.size(),0.0);

  
  double scatterFit=10000.0;

  /**********************************************************

    Test the scattering value of all possible mixture combinations.

   *********************************************************/

  // first read in the mixture list

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
      std::cout<<line<<"\n";
      while(ss.eof()==false){
	ss>>index;
	std::cout<<"index ? "<<index<<"\n";
	mixtureSet.push_back(index);
      }
      mixtureList.push_back(mixtureSet);
      std::cout<<mixtureList.size()<<"\n";
    }
  }else{
      std::cout<<"failed to open mixture file\n";
  }
  
  //now calculate the scattering curve using the implied mixtures
  
  for(int i=0;i<mixtureList.size();i++){
    double scatterFitTemp = ed.fitToScatteringMultiple(molDists,solDists,solMolDists,molSize,noSol,mixtureList[i]);
    std::cout<<"percentage species 1,2,3...: ";
    for(int j=0;j<mixtureList[i].size();j++){
      std::cout<<mixtureList[i][j]<<" ";
    }
    std::cout<<" gives scatter "<<scatterFitTemp<<"\n";
    if(scatterFitTemp<scatterFit){
      scatterFit = scatterFitTemp;
      percentageCombinations = mixtureList[i];
    }
    std::stringstream ss;
    ss<<0;
    const char* str = ss.str().c_str(); 
    char outputLoc[100];
    char absWrFile[100];
    char wrFile[100];
    char scatterFile[100];
    char molFile[100];
    strcpy(outputLoc,argv[13]);
    strcat(outputLoc,"/");
    strcpy(absWrFile,outputLoc);
    strcpy(wrFile,outputLoc);
    strcpy(scatterFile,outputLoc);
    strcpy(molFile,outputLoc);
    strcat(absWrFile,"absWrFP");
    strcat(absWrFile,str);
    strcat(absWrFile,".dat");
    strcat(wrFile,"wrFP");
    strcat(wrFile,str);
    strcat(wrFile,".dat");
    strcat(scatterFile,"scatter");
    strcat(scatterFile,str);
    strcat(scatterFile,".dat");
    strcat(molFile,"mol");
    strcat(molFile,str);
    strcat(molFile,".dat");
    std::cout<<"file ? "<<scatterFile<<"\n";
    std::vector<std::vector<point> > crds =mol[0].getCoordinates();
    localWrithe lw;
    lw.DIDownSampleAbsWrite(crds,absWrFile);
    lw.DIDownSampleWrite(crds,wrFile);
    ed.writeRawMolScatteringToFileMultiple(molDists,solDists,solMolDists,molSize,noSol,mixtureList[i],scatterFile);
    mol[0].writeMoleculeToFile(molFile);
  }
  std::cout<<"done initial scatter\n";

  /***************************************************************

   apply penalties which are "un protein like". Currently we are using

     i) a very strict overlap penalty which exponetiallp penalises non local sections coming close than 4 A.
     ii) A distance constraint measure, which is only active if the user inputs a set of distance consrtrainst like contact predictions.
     iii) A hydrophic covering measure which is very much work in progress but the aim is to make sure there isn't excessive hydrophbic expose 
      to be implemented soon iv) A writhe penalty to ensure the moelule doesn't become too disentangled.
  
   **************************************************************/
  
  
  // Penalise excessive overlap

  double overlapPenalty = 0.0;

  for(int i=0;i<mol.size();i++){
    overlapPenalty = overlapPenalty + getOverlapPenalty(closestApproachDist,overlapDistSet[i]);
  }

  scatterFit =scatterFit + overlapPenalty;

  // add any additional constrainst to the overall fitness based on distance constraints (should rename as its not necssarily Lennard-Jones)

  double contactPredPen=0.0;
  for(int i=0;i<mol.size();i++){
    contactPredPen =contactPredPen + mol[i].getLennardJonesContact();
  }

  scatterFit = scatterFit +  contactPredPen;
  
  double hydrationPenalisation=0.0;
  
  // penalise any excessive hydraion packing (too many solvents close to a hydrophobic amino acid.

  for(int i=0;i<mol.size();i++){
    hydrationPenalisation =hydrationPenalisation + getHydrophobicPackingPenalty(hydrophobicPackingMeasures[i]);
  }

  std::cout<<"initial  hydration penalisation "<<hydrationPenalisation<<"\n";
  scatterFit = scatterFit +  hydrationPenalisation;


  // calculate writhe lists

  std::vector<double> originalWrithes;
  
  for(int i=0;i<mol.size();i++){
   localWrithe lw;
   std::vector<std::vector<point> > crds =mol[i].getCoordinates();
   std::vector<std::pair<std::pair<int,int>,double> > wrFingerPrint =lw.DIDownSampleAbs(crds);
   std::cout<<"wr is "<<wrFingerPrint[wrFingerPrint.size()-1].second<<"\n";
   originalWrithes.push_back(wrFingerPrint[wrFingerPrint.size()-1].second);
  }

  int saveIndex=0;
  
  /******************************************************************
  
                      Fit to scattering data
  
  ******************************************************************/


  /****************************************************************************
    
    Main algorithm 
  
  ***************************************************************************/

  //set up loop parameters
  int k=0;
  
  /* the vector noSections tells us how many subsections are in each moelcule
     e.g. for a monomer/dimer mixture noSections[0]=1,noSections[1]=2.
   */
  
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

    // loop over the molecules

    for(int l=0;l<mol.size();l++){
      int netIndex=0;

      //loop over the sections of the given molecule
      
      for(int i=1;i<=noSections[l];i++){
	
	// net index tells us how far we are through the whole moelcule
	if(i>1){
	  netIndex=netIndex+mol[l].getSubsecSize(i-1);
	}
	// Now loop over the secondary structures of the given unit or section

	
	// make a copy of the molecule (l)  we are looking at we are going to try all chanegs to this and see if changing it improves the best prediction overall.
	ktlMolecule molBest = mol[l];
	/* we will need to store the "best distances" if we find an improved fit changing moelcule l
	   In addition we store best verisons of things like overla, hyrdation fit, best percentage combination e.t.c.
         */
	std::vector<std::vector<double> > molDistsBest = molDists;
	std::vector<std::vector<double> > solDistsBest = solDists;
	std::vector<std::vector<double> > solMolDistsBest = solMolDists;
	std::vector<int> noSolBest = noSol;
	double overlapPenaltyBest=overlapPenalty;
	double hydrationPenalisationBest =hydrationPenalisation;
	int noDistBinsBest;
	std::vector<double> maxDistListBest = maxDistList;
	std::vector<double> hydrophobicPackingMeasuresBest = hydrophobicPackingMeasures;
	std::vector<double> percentageCombinationsBest = percentageCombinations;

	/*********

		  Here wer loop over a subsection of moelcule l (it could be a dimer)

	 *******/
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
	    
	    // this (checkCalphas) checks if there haven't been any rouge sections created (some occasional flaws in the procedure which are to be ironed out
	    bool cacaDist= molCopy.checkCalphas(i);
	    
	    // we are going to change the structure so we will need new distance sets
	    std::vector<std::vector<double> > molDistsTmp = molDists;
	    std::vector<std::vector<double> > solDistsTmp = solDists;
	    std::vector<std::vector<double> > solMolDistsTmp = solMolDists;
	    std::vector<int> noSolTmp = noSol;


	    if(cacaDist==false){
	      // if the new structure is okay we calculate its
	      // generate a new hydration shell
	      std::vector<double> solMolDistsCopy;
	      std::vector<double> solDistsCopy;
	      hydrationShellMinimal hydrationShellCopy(molCopy,Rin,Rout,RShell,ntrivs,helixRatio,solventsPerLink,closestApproachDist,rmin,rmax,lmin);
	      // generate the hydration shell
	      hydrationShellCopy.tubeParamList();
	      hydrationShellCopy.constructInitialState();
	      hydrationShellCopy.getAllHelices();
	      //hydrationShell.allOverlap();
	      hydrationShellCopy.solventMoleculeDistances(solMolDistsCopy,solDistsCopy);
 
	      std::vector<double> molDistsCopy;
	      std::vector<double> meanHydrophilicDistsTmp;
	      //
	      int molSizeTmp,noSolCopy;
	      double maxDistCopy;
	      
	      // check for any overlaps

	      std::vector<double> overlapDistsNew= molCopy.checkOverlapWithRad(closestApproachDist);
	      double newHydrationPackingValue;

	      // generate the new distances
	      generateAllDistances(molCopy,hydrationShellCopy,molDistsCopy,solDistsCopy,solMolDistsCopy,meanHydrophilicDistsTmp,molSizeTmp,noSolCopy,maxDistCopy,newHydrationPackingValue);

	      std::vector<double> maxDistListTmp = maxDistList;
	      std::vector<double> hydrophobicPackingMeasuresTmp = hydrophobicPackingMeasures;
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
		// update hyrophobic packing
		hydrophobicPackingMeasuresTmp[l] = newHydrationPackingValue;
		molDistsTmp[l] = molDistsCopy;
		solDistsTmp[l] = solDistsCopy;
		solMolDistsTmp[l] = solMolDistsCopy;
		noSolTmp[l] = noSolCopy;
		double newScatterFit=10000.0;
		std::vector<double> percentageCombinationsTmp(mol.size(),0.0);
		for(int iv=0;iv<mixtureList.size();iv++){
		  double scatterFitTemp = ed.fitToScatteringMultiple(molDistsTmp,solDistsTmp,solMolDistsTmp,molSize,noSolTmp,mixtureList[iv]);
		  /*std::cout<<"percentage species 1,2,3...: ";
		  for(int jv=0;jv<mixtureList[iv].size();jv++){
		    std::cout<<mixtureList[iv][jv]<<" ";
		  }
		  std::cout<<scatterFitTemp<<"\n";*/
		  if(scatterFitTemp<newScatterFit){
		    newScatterFit = scatterFitTemp;
		    percentageCombinationsTmp = mixtureList[iv];
		  }
		}
		// apply the overlap penalties, only changing for one moelecule hence the if m==l
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
		//std::cout<<" test ? "<<l<<" "<<newScatterFit<<" "<<overlapPenaltyNew<<"\n";
		newScatterFit =newScatterFit +  overlapPenaltyNew;
	      // finally add contract prediction value
		newScatterFit =newScatterFit + contactPredPenNew;
		double hydrationPenalisationNew=0.0;
		
		// penalise any excessive hydration packing (too many solvents close to a hydrophobic amino acid.
		for(int iv=0;iv<mol.size();iv++){
		  hydrationPenalisationNew = hydrationPenalisationNew + getHydrophobicPackingPenalty(hydrophobicPackingMeasuresTmp[l]);
		}//std::cout<<"updated hydration penalisation "<<hydrationPenalisationNew<<"\n";
		newScatterFit = newScatterFit +  hydrationPenalisationNew;

		//
		localWrithe lw;
		std::vector<std::vector<point> > crds =molCopy.getCoordinates();
		std::stringstream ss;
		saveIndex++;
		ss<<saveIndex;
		const char* str = ss.str().c_str(); 
		char outputLoc[100];
		char absWrFile[100];
		char wrFile[100];
		char scatterFile[100];
		char molFile[100];
		strcpy(outputLoc,argv[13]);
		strcat(outputLoc,"/");
		strcpy(absWrFile,outputLoc);
		strcpy(wrFile,outputLoc);
		strcpy(scatterFile,outputLoc);
		strcpy(molFile,outputLoc);
		strcat(absWrFile,"absWrFP");
		strcat(absWrFile,str);
		strcat(absWrFile,".dat");
		strcat(wrFile,"wrFP");
		strcat(wrFile,str);
		strcat(wrFile,".dat");
		strcat(scatterFile,"scatter");
		strcat(scatterFile,str);
		strcat(scatterFile,".dat");
		strcat(molFile,"mol");
		strcat(molFile,str);
		strcat(molFile,".dat");
		lw.DIDownSampleAbsWrite(crds,absWrFile);
	        lw.DIDownSampleWrite(crds,wrFile);
		ed.writeRawMolScatteringToFileMultiple(molDistsTmp,solDistsTmp,solMolDistsTmp,molSize,noSolTmp,percentageCombinationsTmp,scatterFile);
		molCopy.writeMoleculeToFile(molFile);
		// check if we have improved overall, if so update the "best" fit, note that this best is based only on changeing this current molecule.
		double uProb = distributionR(generator1);
		if(checkTransition(newScatterFit,scatterFitBest,uProb,k,noScatterFitSteps)){
		  scatterFitBest  = newScatterFit;       
		  molBest=molCopy;
		  molDistsBest = molDistsTmp;
		  solDistsBest = solDistsTmp;
		  solMolDistsBest = solMolDistsTmp;
		  noSolBest = noSolTmp;
		  overlapPenaltyBest=overlapPenaltyNew;
		  hydrationPenalisationBest = hydrationPenalisationNew;
		  noDistBinsBest=newNoDistBins;
		  maxDistListBest = maxDistListTmp;
		  percentageCombinationsBest = percentageCombinationsTmp;
		  hydrophobicPackingMeasuresBest = hydrophobicPackingMeasuresTmp;
		}
	      }
	    }
	  }
	}
	// after we have found what is best for moelcule l we see if we should update the whole ensemble
	double uProb = distributionR(generator1);
	if(checkTransition(scatterFitBest,scatterFit,uProb,k,noScatterFitSteps)){
	  scatterFit  = scatterFitBest;       
	  mol[l]=molBest;
	  molDists = molDistsBest;
	  solDists = solDistsBest;
	  solMolDists = solMolDistsBest;
	  noSol = noSolBest;
	  overlapPenalty=overlapPenaltyBest;
	  hydrationPenalisation=hydrationPenalisationBest;
	  noDistBins=noDistBinsBest;
	  maxDistList = maxDistListBest;
	  percentageCombinations = percentageCombinationsBest;
	  hydrophobicPackingMeasures=hydrophobicPackingMeasuresBest;
	}
      }
    }
    k++;
    std::cout<<k<<" "<<scatterFit<<" "<<overlapPenalty<<" "<<hydrationPenalisation<<" "<<percentageCombinations[0]<<" "<<percentageCombinations[1]<<"\n";
  }
}
      
