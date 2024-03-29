#include "moleculeFitAndState.h"


moleculeFitAndState::moleculeFitAndState(std::vector<ktlMolecule> &molin,double &RinSv,double &RoutSv,double &RShellSv,int &ntrivsSv,double &closestApproachDistIn,int &solventsPerLinkIn,double &rminIn,double &rmaxIn,double &lminIn){
  //define the number of structures
  mol = molin;
  molDists.resize(mol.size());
  solDists.resize(mol.size());
  solMolDists.resize(mol.size());
  molSize.resize(mol.size());
  noSol.resize(mol.size());
  maxDistMol.resize(mol.size());
  maxDistSol.resize(mol.size());
  contactPredPen.resize(mol.size());
  writhePenalty=0.0;
  // set the fixed fitting parameters
  Rin= RinSv;
  Rout=RoutSv;
  RShell = RShellSv;
  ntrivs=ntrivsSv;
  closestApproachDist = closestApproachDistIn;
  rmin = rminIn;rmax = rmaxIn;solventsPerLink=solventsPerLinkIn;
  lmin = lminIn;
  for(int i=0;i<mol.size();i++){
   writhe lw;
   std::vector<std::vector<point> > crds =mol[i].getCoordinates();
   double acn = lw.DIDownSampleAbs(crds);
  //  std::cout<<"initial abs writhe "<< acn <<"\n";
   originalWrithes.push_back(acn);
  }
  currWrithes =originalWrithes;
  for(int i=0;i<mol.size();i++){
    //calculate distances
    std::vector<double> overlapDists= mol[i].checkOverlapWithRad(closestApproachDist);
    molDists[i] = mol[i].getDistSet();
    //get number of amino acids
    molSize[i]  = mol[i].getNoAminos();
    //sort the distances from largest to smallest for binning.
    std::sort(molDists[i].begin(),molDists[i].end());
    maxDistMol[i] = molDists[i][molDists[i].size()-1];
    // calculate any overlap distances
    overlapDistSet.push_back(overlapDists);
  }
}

std::vector<ktlMolecule> moleculeFitAndState::getMolecule(){
  return mol;
}

void  moleculeFitAndState::updateMolecule(std::vector<ktlMolecule> &molNew){
  mol=molNew;
}

// version where a single changed section has been made and we update the distances

void moleculeFitAndState::calculateMoleculeDistances(ktlMolecule &molNew,int &i){
  std::vector<double> overlapDists= molNew.checkOverlapWithRad(closestApproachDist);
  molDists[i] = molNew.getDistSet();
  //get number of amino acids
  molSize[i]  = molNew.getNoAminos();
  //sort the distances from largest to smallest for binning.
  std::sort(molDists[i].begin(),molDists[i].end());
  maxDistMol[i] = molDists[i][molDists[i].size()-1];
  // calculate any overlap distances
  overlapDistSet[i]=overlapDists;
}


void moleculeFitAndState::calcuateHydrationDistances(hydrationShellMinimal &hs,int &i){
  // we are filling the sol mol distance and sol distance vectors, to prevent double counting we empty it first
  solMolDists[i].clear();
  solDists[i].clear();
  hs.solventMoleculeDistances(solMolDists[i],solDists[i]);
  std::vector<std::vector<point> > solptsin = hs.returnFlatSolList();
  // calculate solvent-solvent distances
  std::sort(solDists[i].begin(),solDists[i].end());
  maxDistSol[i] = solDists[i][solDists[i].size()-1];
  //std::cout<<i<<" "<<solDists[i].size()<<" "<<maxDist<<"\n";
  /* double solMax=5.5;
  hydroPhobicPacking =0.0;
  std::vector<double> meanHydrophilicDists = molin.getHydrophobicDistance(solptsin,solMax);
  for(int i=0;i<meanHydrophilicDists.size();i++){
    hydroPhobicPacking  =  hydroPhobicPacking  + meanHydrophilicDists[i];
  }*/
  //std::cout<<"how much hydro "<<totalHydroOverlap<<"\n";
  noSol[i]=0;
  for(int l=0;l<solptsin.size();l++){
    noSol[i] = noSol[i] + solptsin[l].size();
  }
  //std::cout<<"no amino acis "<<molSizein<<" no solvents "<<noSolin<<"\n";
  std::sort(solMolDists[i].begin(),solMolDists[i].end());
}



// Here we loop over the mixture list (i.e. we do  0.5% monomer, then 0.5% dimer)
// then for each mixture we loop over all hyrdation shell density parametrs, the higher this is, well, the more solvent molcules in the layer
// finall we move over all molecules themselves

double moleculeFitAndState::calculateScattering(experimentalData &ed,double &kmin,double &kmax,std::vector<double> &mixtureVals){
  int noDistBins =ed.setPhases(maxDist,kmin,kmax);
  return  ed.fitToScatteringMultiple(molDists,solDists,solMolDists,molSize,noSol,mixtureVals);
}


void  moleculeFitAndState::writeScatteringToFile(experimentalData &ed,double &kmin,double &kmax,const char* filename){
  int noDistBins =ed.setPhases(maxDist,kmin,kmax);
  ed.writeScatteringToFileMultiple(molDists,solDists,solMolDists,molSize,noSol,percentageCombinations,filename);
}



double moleculeFitAndState::getOverlapPenalty(double &closestApproachDist,std::vector<double> &overlapDists){
  double distSumCurr=0.0;
  for(int l=0;l<overlapDists.size();l++){
    // Uncheck me to see what overalps we get...
    //std::cout<<l<<" "<<overlapDists[l]<<"\n";
    distSumCurr = distSumCurr + std::exp(std::abs(closestApproachDist-overlapDists[l]))-1.0;
  }
  if(overlapDists.size()>0){
    distSumCurr =0.005*distSumCurr/overlapDists.size();
  }
  return distSumCurr;
}

double moleculeFitAndState::applyOverlapPenalty(){
   double overlapPenalty = 0.0;
   for(int i=0;i<overlapDistSet.size();i++){
     // calculate any overlap distances
     overlapPenalty = overlapPenalty + getOverlapPenalty(closestApproachDist,overlapDistSet[i]);
  }
   return overlapPenalty; 
}

double moleculeFitAndState::applyDistanceConstraints(){
 double contactPredPenTotal=0.0;
  for(int i=0;i<mol.size();i++){
    contactPredPen[i] = mol[i].getLennardJonesContact();
    contactPredPenTotal=contactPredPenTotal+contactPredPen[i];
  }
  return contactPredPenTotal;
}

double moleculeFitAndState::applyDistanceConstraints(ktlMolecule &molNew,int &im){
 double contactPredPenTotal=0.0;
  for(int i=0;i<contactPredPen.size();i++){
    if(i==im){
      contactPredPen[i] = molNew.getLennardJonesContact();
      contactPredPenTotal=contactPredPenTotal+contactPredPen[i];
    }else{
      contactPredPenTotal=contactPredPenTotal+contactPredPen[i];
    }
  }
  return contactPredPenTotal;
}

void moleculeFitAndState::alterWrithe(ktlMolecule &molNew,int &i){
  writhe lw;
   std::vector<std::vector<point> > crds =molNew.getCoordinates();
   double acn = lw.DIDownSampleAbs(crds);
   currWrithes[i] = acn;
}

void moleculeFitAndState::alterWritheSet(ktlMolecule &molNew,int &i){
  writhe lw;
   std::vector<std::vector<point> > crds =molNew.getCoordinates();
   std::vector<std::pair<std::pair<int,int>,double> > wrFingerPrint =lw.DIDownSampleAbsFP(crds);
   currWrithes[i] = wrFingerPrint[wrFingerPrint.size()-1].second;
}

// calculate writhe lists
// void moleculeFitAndState::applyWritheConstraint(){
//   writhePenalty=0.0;
//   for(int i=0;i<originalWrithes.size();i++){
//    double newWrithe =currWrithes[i];
//    double secLen = double(mol[i].molSize());
//    double lowerBound = std::pow((secLen/7.5),1.625);
//    writhePenalty=  writhePenalty+1.0/(1.0+std::exp(-(newWrithe-lowerBound)));
//   }
// }



void moleculeFitAndState::applyWritheConstraint(){

  writhePenalty=0.0;

  for(int i=0;i<originalWrithes.size();i++){

   double newWrithe =currWrithes[i];

   double secLen = double(mol[i].noSecSize());

   double lowerBound = std::pow((secLen/7.5),1.625);

   writhePenalty=  writhePenalty+1.0/(1.0+std::exp(20.0*(newWrithe-lowerBound)));

  }

}

double moleculeFitAndState::getFit(){
  return currFit;
}



std::tuple<double, double, double, double> moleculeFitAndState::getOverallFit(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,std::vector<double> &helRatList,double &kmin,double &kmax){
  double scatterAndHydrationConstraint =10000.0;
  int bestHelRatList=0;
  for(int m=0;m<mixtureList.size();m++){
    for(int j=0;j<helRatList.size();j++){
      for(int i=0;i<mol.size();i++){
	//generate the hydration shell
	hydrationShellMinimal hydrationShellTmp(mol[i],Rin,Rout,RShell,ntrivs,helRatList[j],solventsPerLink,closestApproachDist,rmin,rmax,lmin);
	hydrationShellTmp.generateHydrationLayer();
	//calclate the distances associated with the shell
	calcuateHydrationDistances(hydrationShellTmp,i);
	// get the scattering value
	//std::cout<<m<<" "<<j<<" "<<i<<" "<<mixtureList.size()<<" "<<helRatList.size()<<" "<<mol.size()<<"\n";
      }
      double tempScat=11000.0;
      double maxDistMolecule = *std::max_element(maxDistMol.begin(),maxDistMol.end());
      maxDist = *std::max_element(maxDistSol.begin(),maxDistSol.end());
      //std::cout<<maxDist<<" "<<maxDistMolecule<<"\n";
      if(maxDist<2.0*maxDistMolecule){
	tempScat =  calculateScattering(ed,kmin,kmax,mixtureList[m]);
      }
      //check what is best
      //std::cout<<tempScat<<"\n";
      if(tempScat< scatterAndHydrationConstraint){
	scatterAndHydrationConstraint= tempScat;
	percentageCombinations = mixtureList[m];
	bestHelRatList=j;
      }
    }  
  }

  // regenerate the best one
  for(int i=0;i<mol.size();i++){
    //generate the hydration shell
    hydrationShellMinimal hydrationShellTmp(mol[i],Rin,Rout,RShell,ntrivs,helRatList[bestHelRatList],solventsPerLink,closestApproachDist,rmin,rmax,lmin);
    hydrationShellTmp.generateHydrationLayer();
    //calclate the distances associated with the shell
    calcuateHydrationDistances(hydrationShellTmp,i);
  }
  /***************************************************************

   apply penalties which are "un protein like". Currently we are using

     i) a very strict overlap penalty which exponetiallp penalises non local sections coming close than 4 A.
     ii) A distance constraint measure, which is only active if the user inputs a set of distance consrtrainst like contact predictions.
     iii) A writhe penalty to ensure the moelule doesn't become too disentangled.
  
  **************************************************************/
  double overlapPenalty = applyOverlapPenalty();
  // std::cout<<"overlap penalty "<<overlapPenalty<<"\n";
  double distanceConstraints = applyDistanceConstraints();
  // std::cout<<"distance constraints "<<distanceConstraints<<"\n";
  applyWritheConstraint();
  // std::cout<<"writhe penalty "<<writhePenalty<<"\n";
  // std::cout<<" scattering  "<<scatterAndHydrationConstraint<<"\n";
  currFit = scatterAndHydrationConstraint +overlapPenalty +distanceConstraints + writhePenalty;
  return std::make_tuple(currFit, overlapPenalty, writhePenalty, scatterAndHydrationConstraint);
}


std::tuple<double, double, double, double> moleculeFitAndState::getOverallFit(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,std::vector<double> &helRatList,ktlMolecule &molNew,double &kmin,double &kmax,int &i){
  // update the molecule distances for molecule i;
  calculateMoleculeDistances(molNew,i);
  int bestHelRatList=0;
  // now update the hydration shell
  double scatterAndHydrationConstraint =10000.0;
  for(int m=0;m<mixtureList.size();m++){
    for(int j=0;j<helRatList.size();j++){
	//generate the hydration shell
	hydrationShellMinimal hydrationShellTmp(molNew,Rin,Rout,RShell,ntrivs,helRatList[j],solventsPerLink,closestApproachDist,rmin,rmax,lmin);
	hydrationShellTmp.generateHydrationLayer();
	//calclate the distances associated with the shell
	calcuateHydrationDistances(hydrationShellTmp,i);
	// get the scattering value
	//std::cout<<m<<" "<<j<<" "<<i<<" "<<mixtureList.size()<<" "<<helRatList.size()<<" "<<mol.size()<<"\n";
	double tempScat=11000.0;
	double maxDistMolecule = *std::max_element(maxDistMol.begin(),maxDistMol.end());
	maxDist = *std::max_element(maxDistSol.begin(),maxDistSol.end());
	//std::cout<<maxDist<<" "<<maxDistMolecule<<"\n";
	if(maxDist<2.0*maxDistMolecule){
	  tempScat =  calculateScattering(ed,kmin,kmax,mixtureList[m]);
	}
	//check what is best
	//std::cout<<tempScat<<"\n";
	if(tempScat< scatterAndHydrationConstraint){
	  scatterAndHydrationConstraint= tempScat;
	  percentageCombinations = mixtureList[m];
	  bestHelRatList=j;
	}
    }
  }
  // regenerate the best one
  hydrationShellMinimal hydrationShellTmp(molNew,Rin,Rout,RShell,ntrivs,helRatList[bestHelRatList],solventsPerLink,closestApproachDist,rmin,rmax,lmin);
  hydrationShellTmp.generateHydrationLayer();
  //calclate the distances associated with the shell
  calcuateHydrationDistances(hydrationShellTmp,i);
  // apply penalties
   double overlapPenalty = applyOverlapPenalty();
  // std::cout<<"Overlap Penalty "<<overlapPenalty<<"\n";
   double distanceConstraints = applyDistanceConstraints(molNew,i);
   //std::cout<<"pen3 "<<distanceConstraints<<"\n";
  alterWrithe(molNew,i);
  applyWritheConstraint();
  // std::cout<<" writhe penalty  "<<writhePenalty<<"\n";
  // std::cout<<" scattering  "<<scatterAndHydrationConstraint<<"\n";
  currFit = scatterAndHydrationConstraint +overlapPenalty +distanceConstraints + writhePenalty;
  return std::make_tuple(currFit, overlapPenalty, writhePenalty, scatterAndHydrationConstraint);
}
