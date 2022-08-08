
#include "hydrationShellRandom.h"



hydrationShellMinimal::hydrationShellMinimal(ktlMolecule &molIn,double RinInp,double RoutInp,double &RhyIn,int ntrivsIn,double helixRatioIn,int solventsPerLinkIn,double &mutualDistCutOffIn,double &rmin,double &rmax,double &lmin){
  mol = molIn;
  molsize = mol.molSize();
  nSize = mol.noSecSize();
  direcList.resize(molsize);
  frameTan.resize(molsize);
  frameNorm.resize(molsize);
  frameBinorm.resize(molsize);
  halfLengthList.resize(molsize);
  nameLst.resize(molsize);
  lengthSec.resize(molsize);
  midPointList.resize(molsize);
  Pi = 3.1415926535897932384626;
  edgesolset1.resize(3);
  edgesolset2.resize(3);
  edgesolset3.resize(3);
  edgesolset4.resize(3);
  middlesol.resize(3);
  Rin = RinInp;Rout = RoutInp;Rhy =RhyIn;
  allSegments.resize(molsize);
  allTruthTables.resize(molsize);
  solventsPerLink = solventsPerLinkIn;
  helixRatio = helixRatioIn;
  ntrivs = ntrivsIn;
  helixpts.resize(molsize);
  cstepList.resize(molsize);
  selfScatteringListSol.resize(molsize);
  selfScatteringListMol.resize(molsize);
  selfScatteringListSolMol.resize(molsize);
  mutualScatteringListSol.resize(molsize);
  mutualScatteringListMol.resize(molsize);
  mutualScatteringListSolMol.resize(molsize);
  for(int i=0;i<molsize;i++){
    mutualScatteringListSol[i].resize(molsize);
    mutualScatteringListMol[i].resize(molsize);
    mutualScatteringListSolMol[i].resize(molsize);
    std::vector<std::pair<double,int> > distPosList;
    if(i<molsize-1){
      for(int j=i+1;j<molsize;j++){
      std::pair<double,int> p;
      p.first = 1000000.0;
      p.second = j;
    	distPosList.push_back(p);
    }
      std::pair<std::vector<std::pair<double,int> >,int> minpr;
      minpr.first = distPosList;
      minpr.second =0;
      minimumDistances.push_back(minpr);
    }
  }
  hydrationChangeList.resize(molsize);
  //declare the empty distances list
  maxDistChange =0.0;
  std::pair<double,double> maxPair=mol.getMaxPossibleLength();
  maxdist=maxPair.first;
  maxScatDist = maxPair.second;
  mutualDistCutOff =mutualDistCutOffIn;
  isTooClose= false;
}


void hydrationShellMinimal::resetMol(ktlMolecule &molIn){
  mol = molIn;
}


point hydrationShellMinimal::getCentrePoint(int index,int subIndex){
  return mol.getCoordinate(index,subIndex);
}

int hydrationShellMinimal::getNoKvals(){
  return nscat;
}

double hydrationShellMinimal::getMaxDist(){
  return maxdist;
}  

void hydrationShellMinimal::getPointAndMidlengthHelix(int i,int &hIndex){
  double kappa = mol.getCurvatureJoined(i);
  double tau = mol.getTorsionJoined(i);
  int noMols = mol.getUnitNo(i);
  double L = mol.getAlbeadJoined(i)*(noMols-1);
  point tan = mol.getTangent(i,0);
  point norm = mol.getNormal(i,0);
  point binorm = mol.getBinormal(i,0);
  point startCord = mol.getCoordinate(i,0);
  double rad = kappa/(kappa*kappa + tau*tau);
  //start point
  double sv = 0.0;
  ph.updateYvecGeneral(kappa,tau,tan,norm,binorm,startCord,sv);
  point startPt = ph.getCoord() + ph.getNormal()*rad;
  //mid point 
  sv =0.5*L;
  ph.updateYvecGeneral(kappa,tau,tan,norm,binorm,startCord,sv);  
  point midPt = ph.getCoord() + ph.getNormal()*rad;
  // at the mid point we also calculate the frame and the mid point direction
  point midtan = ph.getTangent();
  point midNorm = ph.getNormal();
  point midBinorm = ph.getBinormal();
  point direc = midBinorm*kappa + midtan*tau;
  midPointList[hIndex]=midPt;
  frameTan[hIndex] = midtan;
  frameNorm[hIndex] = midNorm;
  frameBinorm[hIndex] = midBinorm;
  //end point
  sv = L;
  ph.updateYvecGeneral(kappa,tau,tan,norm,binorm,startCord,sv);  
  point endPoint = ph.getCoord() + ph.getNormal()*rad;
  if(L>0.0000001){
    double endDist = startPt.eDist(endPoint);
    double direcWeight = 2.0*Pi*tau/std::pow(kappa*kappa+tau*tau,1.5);
    direc = direc*direcWeight;
    direc.normalise();
    direcList[hIndex]=direc;
    halfLengthList[hIndex]=endDist*0.5;
  }else{
    point p(0.0,0.0,0.0);
    direcList[hIndex]=p;  
    halfLengthList[hIndex]=0.0;
  }
}

void hydrationShellMinimal::getPointAndMidlengthMulti(int i,int &hIndex){
  std::vector<point> coords  = mol.getCoordinatesSection(i);int sz =coords.size();
  std::vector<point> lines;point mean(0.0,0.0,0.0);
  for(int j=0;j<sz;j++){
    mean = mean + coords[j];
  }
  mean = mean*(1.0/double(sz));
  point endDif = mean - coords[sz-1];
  for(int j=0;j<sz;j++){
    lines.push_back(mean - coords[j]);
    if(lines[j].dotprod(endDif)<0.0){
      lines[j] = lines[j]*(-1.0);
    }
  }
  point direc(0.0,0.0,0.0);
  for(int j=0;j<sz;j++){
    direc = direc + lines[j];
  }
  double hlfDist= 0.5*(mean.eDist(coords[0])+mean.eDist(coords[sz-1]));
  direc.normalise();
  midPointList[hIndex]=mean;
  frameTan[hIndex] = direc;
  point norm;
  point vert(0.0,0.0,1.0);
  if(std::abs(direc.dotprod(vert))<0.99999999){
    norm = vert.cross(direc);
    norm.normalise();
  }else{
    point vert2(0.0,1.0,0.0);
    norm= vert2.cross(direc);
    norm.normalise();
  }
  frameNorm[hIndex] = norm;
  frameBinorm[hIndex] = direc.cross(norm);
  direcList[hIndex]=direc;
  halfLengthList[hIndex]=hlfDist;
}


void hydrationShellMinimal::getPointAndMidlengthStraight(int &sec,int &part,int &hIndex,std::string soe){
  point startPt,endPt,midPoint;
  if(soe == "start"){
    if(sec==0){
      midPoint = mol.getCoordinate(sec,part);
      endPt = mol.getCoordinate(sec,part+1);
      startPt = midPoint -(endPt-midPoint);
    }else{    
      startPt = mol.getCoordinate(sec,-1);
      midPoint = mol.getCoordinate(sec,part);
      endPt = mol.getCoordinate(sec,part+1);
    }
  }else if(soe == "end"){
    if(sec==nSize-1){
      startPt = mol.getCoordinate(sec,part-1);
      midPoint = mol.getCoordinate(sec,part);
      endPt = midPoint+(midPoint-startPt);
    }else{
      startPt = mol.getCoordinate(sec,part-1);
      midPoint = mol.getCoordinate(sec,part);
      endPt = mol.getCoordinate(sec+1,0);
    }
  }else{
    startPt = mol.getCoordinate(sec,part-1);
    midPoint = mol.getCoordinate(sec,part);
    endPt = mol.getCoordinate(sec,part+1);
  }
  double endDist = 0.5*startPt.eDist(endPt);
  point direc = endPt-startPt;
  direc.normalise();
  direcList[hIndex]=direc;
  frameTan[hIndex] = direc;
  double normt;
  if(sec==0||sec == nSize-1){
    point vert(0.0,0.0,1.0);
    if(std::abs(direc.dotprod(vert))<0.9999999){
      point N1 = vert.cross(direc);
      N1.normalise();
      frameNorm[hIndex] = N1;
    }else{
      point vert2(0.0,1.0,0.0);
      point N1 = vert2.cross(direc);
      N1.normalise();
      frameNorm[hIndex] = N1;
    }
  }else{
    point bmina = endPt-startPt;
    normt = (midPoint-startPt).dotprod(bmina);
    normt = normt/(bmina.dotprod(bmina));
    point frnm =  startPt + (endPt-startPt)*normt-midPoint;
    frnm.normalise();
    frameNorm[hIndex] =frnm;
  }
  frameNorm[hIndex].normalise();
  frameBinorm[hIndex] = frameTan[hIndex].cross(frameNorm[hIndex]);
  midPointList[hIndex]= startPt;   
  halfLengthList[hIndex]=endDist*0.5;
}


void hydrationShellMinimal::tubeParamList(){
  int currHIndex=0;
  for(int i=0;i<nSize;i++){
    if(mol.getType(i)=="Helix" || mol.getType(i)=="Strand"){
      getPointAndMidlengthMulti(i,currHIndex);
      nameLst[currHIndex] =mol.getType(i);
      lengthSec[currHIndex] = mol.getUnitNo(i);
      currHIndex++;
    }else{
      // loop over the whole section
      for(int j=0;j<mol.getUnitNo(i);j++){
        if(j==0){
	  getPointAndMidlengthStraight(i,j,currHIndex,"start");
	}else if(j==mol.getUnitNo(i)-1){
	  getPointAndMidlengthStraight(i,j,currHIndex,"end");
	}else{
	  getPointAndMidlengthStraight(i,j,currHIndex,"middle");
	}
        nameLst[currHIndex] = "loop";
        lengthSec[currHIndex] = 1;
	currHIndex++;
      }
    }
  }
}

void hydrationShellMinimal::updateTubeParamList(int &index,int &startHIndex){
  // need to calculate what currHindex is
  int currHIndex =startHIndex;
  for(int i=index-1;i<nSize;i++){
    if(mol.getType(i)=="Helix" || mol.getType(i)=="Strand"){
      getPointAndMidlengthMulti(i,currHIndex);
      currHIndex++;
    }else{
      // loop over the whole section
      for(int j=0;j<mol.getUnitNo(i);j++){
	if(j==0){
	  getPointAndMidlengthStraight(i,j,currHIndex,"start");
	}else if(j==mol.getUnitNo(i)-1){
	  getPointAndMidlengthStraight(i,j,currHIndex,"end");
	}else{
	  getPointAndMidlengthStraight(i,j,currHIndex,"middle");
	}
	currHIndex++;
      }
    }
  }
}

void hydrationShellMinimal::getAllHelices(){
  int currHIndex=0;
  for(int i=0;i<nSize;i++){
    if(mol.getType(i)=="Helix" || mol.getType(i)=="Strand"){
      helixpts[currHIndex]=mol.getCoordinatesSection(i);
      currHIndex++;
    }else{
      // loop over the whole section
      for(int j=0;j<mol.getUnitNo(i);j++){
	      std::vector<point> helpt;
      	helpt.push_back(mol.getCoordinate(i,j));
        helixpts[currHIndex] = helpt;
      	currHIndex++;
      }
    }
  }
};


void hydrationShellMinimal::updateHelices(int &changeIndex){
  int currHIndex =0;
  for(int l=0;l<changeIndex-1;l++){
    if(mol.getType(l)=="Helix" || mol.getType(l)=="Strand"){
      currHIndex = currHIndex+1;
    }else{
      currHIndex = currHIndex+mol.getUnitNo(l);
    }
  }
  for(int i=changeIndex-1;i<nSize;i++){
    if(mol.getType(i)=="Helix" || mol.getType(i)=="Strand"){
      helixpts[currHIndex]=mol.getCoordinatesSection(i);
      currHIndex++;
    }else{
      // loop over the whole section
      for(int j=0;j<mol.getUnitNo(i);j++){
	std::vector<point> helpt;
	helpt.push_back(mol.getCoordinate(i,j));
        helixpts[currHIndex] = helpt;
	currHIndex++;
      }
    }
  }
}



point hydrationShellMinimal::getDirec(int i){
  return direcList[i];
}

point hydrationShellMinimal::getHydTan(int i){
  return frameTan[i];
}

point hydrationShellMinimal::getHydNorm(int i){
  return frameNorm[i];
}

point hydrationShellMinimal::getHydBinorm(int i){
  return frameBinorm[i];
}

double hydrationShellMinimal::hydTubeLength(int i){
  return 2*halfLengthList[i];
}

point hydrationShellMinimal::hydTubeCentre(int i){
  return midPointList[i];
}

int hydrationShellMinimal::getMolSize(){
  return mol.molSize();
}

int hydrationShellMinimal::getNoSections(){
  return mol.noSecSize();
}


point hydrationShellMinimal::getTangent(int secindex,int subIndex){
  return mol.getTangent(secindex,subIndex);
};

point hydrationShellMinimal::getNormal(int secindex,int subIndex){
  return mol.getNormal(secindex,subIndex);
};

point hydrationShellMinimal::getBinormal(int secindex,int subIndex){
  return mol.getBinormal(secindex,subIndex);
};

int hydrationShellMinimal::getNScat(){
  return nscat;
}

/*void hydrationShellMinimal::resetMolecule(){
  mol.resetRandomMolecule();
  }*/

void hydrationShellMinimal::makeInitialSegData(point &cp,point &T,point &N1,double &tm,int index,int &nseg){
  point N12 = T.cross(N1);
  double cstep;
  if(nseg>1){
    cstep = 2*tm/double(nseg-1);
  }else{
    cstep =0.0;
  }
  //find the points along the "tube" axis where the hydrations discs will be centered
  cstepList[index] = cstep;
  std::vector<point> line;
  if(nseg%2==0){
    for(int i=-nseg/2;i <= -1;i++){
      line.push_back(cp + T*(i*cstep + 0.5*cstep));
    }
    for(int i=1;i<= nseg/2;i++){
      line.push_back(cp + T*(i*cstep - 0.5*cstep));
    }
  }else{
    if(nseg>1){
      for(int i=-(nseg-1)/2;i<=(nseg-1)/2;i++){
        line.push_back(cp + T*(i*cstep));
      }
    }else{
        line.push_back(cp);
    }
  }
  std::vector<std::vector<point> > segments;
  std::vector<std::vector<int> > truthTable(line.size(),std::vector<int>(ntrivs,1));
  point p;
  std::vector<int> checkList;
  checkList.push_back(1000000);
  std::vector<std::vector<int> > checkDisc;
  std::vector<std::vector<std::vector<int> > > checkFull;
  for(int i=0;i<line.size();i++){
    std::vector<point> seg(ntrivs);
    for(int j=1;j<=ntrivs;j++){
      p = line[i] - T*0.5*cstep + (N1*std::cos(2*Pi*j/ntrivs + Pi/ntrivs) + N12*std::sin(2*Pi*j/ntrivs + Pi/ntrivs))*Rhy;
      seg[j-1] =p;
      checkDisc.push_back(checkList);
    }
    checkFull.push_back(checkDisc);
    segments.push_back(seg);
  }
  allSegments[index]= segments;
  allTruthTables[index] = truthTable;
  // now add the initial index list
  indexList[index]=checkFull;
}



void hydrationShellMinimal::updateSegData(point &cp,point &T,point &N1,double &tm,int index,int &nseg){
  point N12 = T.cross(N1);
  double cstep;
  if(nseg>1){
    cstep = 2*tm/double(nseg-1);
  }else{
    cstep =0.0;
  }
  /*find the points along the "tube" axis where the hydrations discs will be centered*/
  std::vector<point> line;
  if(nseg%2==0){
    for(int i=-nseg/2;i <= -1;i++){
      line.push_back(cp + T*(i*cstep + 0.5*cstep));
    }
    for(int i=1;i<= nseg/2;i++){
      line.push_back(cp + T*(i*cstep - 0.5*cstep));
    }
  }else{
    if(nseg>1){
      for(int i=-(nseg-1)/2;i<=(nseg-1)/2;i++){
        line.push_back(cp + T*(i*cstep));
      }
    }else{
        line.push_back(cp);
    }
  }
  point p;
  std::vector<std::vector<point> > segments;
  for(int i=0;i<line.size();i++){
    std::vector<point> seg(ntrivs);
    for(int j=1;j<=ntrivs;j++){
        p = line[i] - T*0.5*cstep + (N1*std::cos(2*Pi*j/ntrivs + Pi/ntrivs) + N12*std::sin(2*Pi*j/ntrivs + Pi/ntrivs))*Rhy;
      seg[j-1] =p;
    }
    segments.push_back(seg);
  }
  allSegments[index]= segments;
}


void hydrationShellMinimal::constructInitialState(){
  indexList.resize(molsize);
  for(int i=0;i<molsize;i++){
    point outerNorm= frameNorm[i]*(-1.0);
    if(nameLst[i]=="Helix"||nameLst[i]=="Strand"){
      // helix section
      //int noSegs = std::max(1,int(ceil(helixRatio*noCalphas)));
      int noCalphas = lengthSec[i];
      int noSegs = std::max(1,int(round(helixRatio*noCalphas)));
      makeInitialSegData(midPointList[i],direcList[i],outerNorm,halfLengthList[i],i,noSegs);
    }else{
      // linker
      makeInitialSegData(midPointList[i],direcList[i],outerNorm,halfLengthList[i],i,solventsPerLink);
    }
  }
}


void hydrationShellMinimal::updateState(int &index){
  for(int i=index;i<molsize;i++){
    point outerNorm= frameNorm[i]*(-1.0);
     if(nameLst[i]=="Helix"||nameLst[i]=="Strand"){
      // helix section
       int noCalphas = lengthSec[i];
      int noSegs = std::max(1,int(ceil(helixRatio*noCalphas)));
      //int noSegs = std::max(1,int(round(helixRatio*noCalphas)));
      updateSegData(midPointList[i],direcList[i],outerNorm,halfLengthList[i],i,noSegs);
    }else{
      // linker
      updateSegData(midPointList[i],direcList[i],outerNorm,halfLengthList[i],i,solventsPerLink);
    }
  }
}



std::vector<std::vector<point> > hydrationShellMinimal::getHydrationLayer(int i){
  return allSegments[i];
}

std::vector<std::vector<int> > hydrationShellMinimal::getTruthTable(int i){
  return allTruthTables[i];
}


std::pair<std::vector<std::vector<double> >,std::pair<std::string,std::string> >  hydrationShellMinimal::getMinTvalsFiniteNew(point &cp1,point &cp2,point &T1,point &T2,double &tm1,double &tm2){
  std::pair<std::vector<std::vector<double> >,std::pair<std::string,std::string> > output;
  point cpd = cp1 -cp2;
  point cpd2 = cp2 -cp1;
  double tdot = T1.dotprod(T2);
  double t1 = (cpd.dotprod(T2)*tdot - cpd.dotprod(T1))*(1.0/(1.0-tdot*tdot));
  double t2 = (cpd.dotprod(T1)*tdot*(-1.0) + cpd.dotprod(T2))*(1.0/(1.0-tdot*tdot));
  /*double t1mu2 = (cpd2.dotprod(T2) + tm2)*(1.0/tdot);
  double t1ml2 = (cpd2.dotprod(T2) - tm2)*(1.0/tdot);
  double t2mu1 = (c + tm1)*(1.0/tdot);
  double t2ml1 = (cpd.dotprod(T1) - tm1)*(1.0/tdot);*/
  double t1mu2 = tdot*tm2 -cpd.dotprod(T1);
  double t1ml2 = -tdot*tm2 -cpd.dotprod(T1);
  double t2mu1 = tdot*tm1 -cpd2.dotprod(T2);
  double t2ml1 =  -tdot*tm1 -cpd2.dotprod(T2);
  std::vector<double> edgesolset1(3),edgesolset2(3),edgesolset3(3),edgesolset4(3);
  std::vector<std::vector<double> > sols(4); 
  if(t1< -tm1 || t1 > tm1 || t2< -tm2 || t2>tm2){
    /* closest distance is not within the parpameter range, one of the ends will be involved asw the nearest point*/
    if(t1mu2> tm1){
      edgesolset1[0] = tm1; edgesolset1[1] = tm2;
      point pe1 = cp1 + T1*t1;
      point pe2 = cp2 + T2*tm2;
      edgesolset1[2]= pe1.eDist(pe2);
    }else if(t1mu2< -tm1){
      edgesolset1[0] = -tm1; edgesolset1[1] = tm2;
      point pe1 = cp1 - T1*tm1;
      point pe2 = cp2 + T2*tm2;
      edgesolset1[2]= pe1.eDist(pe2);
    }else{
      edgesolset1[0] = t1mu2; edgesolset1[1] = tm2;
      point pe1 = cp1 + T1*t1mu2;
      point pe2 = cp2 + T2*tm2;
      edgesolset1[2]= pe1.eDist(pe2);
    }
    if(t1ml2 >tm1){
      edgesolset2[0] = tm1; edgesolset2[1] = -tm2;
      point pe1 = cp1 + T1*tm1;
      point pe2 = cp2 - T2*tm2;
      edgesolset2[2]= pe1.eDist(pe2);
    }else if(t1ml2 <-tm1){
      edgesolset2[0] = -tm1; edgesolset2[1] = -tm2;
      point pe1 = cp1 - T1*tm1;
      point pe2 = cp2 - T2*tm2;
      edgesolset2[2]= pe1.eDist(pe2);
    }else{
      edgesolset2[0] = t1ml2; edgesolset2[1] = -tm2;
      point pe1 = cp1 + T1*t1ml2;
      point pe2 = cp2 - T2*tm2;
      edgesolset2[2]= pe1.eDist(pe2);
    }
    if(t2mu1>tm2){
      edgesolset3[0] = tm1; edgesolset3[1] = tm2;
      point pe1 = cp1 + T1*tm1;
      point pe2 = cp2 + T2*tm2;
      edgesolset3[2]= pe1.eDist(pe2);
    }else if(t2mu1<-tm2){
      edgesolset3[0] = tm1; edgesolset3[1] = -tm2;
      point pe1 = cp1 + T1*tm1;
      point pe2 = cp2 - T2*tm2;
      edgesolset3[2]= pe1.eDist(pe2);
    }else{
      edgesolset3[0] = tm1; edgesolset3[1] = t2mu1;
      point pe1 = cp1 + T1*tm1;
      point pe2 = cp2 + T2*t2mu1;
      edgesolset3[2]= pe1.eDist(pe2);
    }
    if(t2ml1> tm2){
      edgesolset4[0] = -tm1; edgesolset4[1] = tm2;
      point pe1 = cp1 - T1*tm1;
      point pe2 = cp2 + T2*tm2;
      edgesolset4[2]= pe1.eDist(pe2);
    }else if(t2ml1 <- tm2){
      edgesolset4[0] = -tm1; edgesolset4[1] = -tm2;
      point pe1 = cp1 - T1*tm1;
      point pe2 = cp2 - T2*tm2;
      edgesolset4[2]= pe1.eDist(pe2);
    }else{
      edgesolset4[0] = -tm1; edgesolset4[1] = t2ml1;
      point pe1 = cp1 - T1*tm1;
      point pe2 = cp2 + T2*t2ml1;
      edgesolset4[2]= pe1.eDist(pe2);
    }
    typeMessage.first = "edge";
    typeMessage.second = "edge";
    sols[0] = edgesolset1;sols[1] = edgesolset2;
    sols[2] = edgesolset3;sols[3] = edgesolset4;
    output.first = sols;
    output.second = typeMessage;
  }else{
    if((t1mu2 >= -tm1 && t1mu2 <= tm1) || (t1ml2 >= -tm1 && t1ml2 <= tm1)){
      middlesol[0] = t1;middlesol[1] = t2;
      sols[0] = middlesol;
       output.first = sols;
      typeMessage.first = "middle";
      typeMessage.second = "1contains2";
    }else{
      middlesol[0] = t1;middlesol[1] = t2;
      sols[0] = middlesol;
      output.first = sols;
      typeMessage.first = "middle";
      typeMessage.second = "2contains1";
      output.second = typeMessage;
    }
  }
  return output;
}

bool hydrationShellMinimal::checkIsIn2(point &cp1,point &cp2,point &T1,point &T2,double &tm1,double &tm2,double &sol1,double &sol2,int inch,int jnch,double &R1,double &R2){
    //get the end points
    point pe1 = cp1 + T1*sol1;
    point pe2 = cp2 + T2*sol2;
    // gett he join vector, from  1 to 2
    point J = pe2 - pe1;
    point Jorig = J;
    J.normalise();
    // Project J into the plane
    // need to allow for the case in which J is parallelto T1
    point Jperp = J - T1*J.dotprod(T1);
    Jperp.normalise();
    point pend1 = pe1 + Jperp*R1;
    point pend2 = pe1 - Jperp*R1;
    //check if these points are in 
    double tf1 = (pend1-cp2).dotprod(T2);
    double tf2 = (pend2-cp2).dotprod(T2);
    point jp1 = cp2 + T2*tf1;
    point jp2 = cp2 + T2*tf2;
    double mindist1,mindist2;
    point Jperp2;
    if(tf1<tm2 && tf1>-tm2){
      // nearest point is potentially inside the tube
      mindist1 = jp1.eDist(pend1)-R2;
    }else{
      // look for a new nearets point on the cap of the tube
      Jperp2= Jorig - T2*Jorig.dotprod(T2);
      Jperp2.normalise();
      // project onto then "end plane" of T2
      double rval1 = (pend1-pe2).dotprod(Jperp2);
      if(std::abs(rval1) > R2){
	      // here the nearest point on the "end plane" is out of the allowed radius, find distance to the nearest edge
      	if(rval1 >0){
           point pt = pe2+Jperp2*R2;
           mindist1 = pend1.eDist(pt);
	      }else{
           point pt = pe2-Jperp2*R2;
           mindist1 = pend1.eDist(pt);
	      }
      }else{
	        //here the nearest poin in the end plane is within its radius
	        point pt = pe2+Jperp2*rval1;
	        mindist1 = pend1.eDist(pt);
      }
    }
    // for the second possibility
    if(tf2<tm2 && tf2>-tm2){
      // nearest point is potentially inside the tube
      mindist2 = jp2.eDist(pend2)-R2;
    }else{
      // project onto then "end plane" of T2
      double rval2 = (pend2-pe2).dotprod(Jperp2);
      if(std::abs(rval2) > R2){
	      // here the nearest point on the "end plane" is out of the allowed radius, find distance to the nearest edge
      	if(rval2 >0){
           point pt = pe2+Jperp2*R2;
           mindist2 = pend2.eDist(pt);
	      }else{
           point pt = pe2-Jperp2*R2;
           mindist2 = pend2.eDist(pt);
	      }
      }else{
	        //here the nearest poin in the end plane is within its radius
	        point pt = pe2+Jperp2*rval2;
	        mindist2 = pend2.eDist(pt);
      }
    }
    // update min distance (if necessary) and cheeck of we have found overlap
    if(mindist1>mindist2){
      if(mindist2<minimumDistances[inch].first[jnch-inch-1].first){
	minimumDistances[inch].first[jnch-inch-1].first =mindist2;  
      }
      if(mindist2<0.0){
        return true;
      }else{
        return false;
      }
    }else{
      if(mindist1<minimumDistances[inch].first[jnch-inch-1].first){
	minimumDistances[inch].first[jnch-inch-1].first =mindist1;  
      }
      if(mindist1<0.0){
        return true;
      }else{
        return false;
      }
    }
}


bool hydrationShellMinimal::checkIsIn1(point &cp1,point &cp2,point &T1,point &T2,double &tm1,double &tm2,double &sol1,double &sol2,int inch,int jnch,double &R1,double &R2){
    //get the end points
    point pe1 = cp1 + T1*sol1;
    point pe2 = cp2 + T2*sol2;
    // gett he join vector, from  1 to 2
    point J = pe1 - pe2;
    point Jorig = J;
    J.normalise();
    // Project J into the plane
    // need to allow for the case in which J is parallelto T2
    point Jperp = J - T2*J.dotprod(T2);
    Jperp.normalise();
    point pend1 = pe2 + Jperp*R2;
    point pend2 = pe2 - Jperp*R2;
    //check if these points are in 
    double tf1 = (pend1-cp1).dotprod(T1);
    double tf2 = (pend2-cp1).dotprod(T1);
    point jp1 = cp1 + T1*tf1;
    point jp2 = cp1 + T1*tf2;
    double mindist1,mindist2;
    point Jperp2;
    if(tf1<tm1 && tf1>-tm1){
      // nearest point is potentially inside the tube
      mindist1 = jp1.eDist(pend1)-R1;
    }else{
      // look for a new nearets point on the cap of the tube
      Jperp2= Jorig - T1*Jorig.dotprod(T1);
      Jperp2.normalise();
      // project onto then "end plane" of T2
      double rval1 = (pend1-pe1).dotprod(Jperp2);
      if(std::abs(rval1) > R1){
	      // here the nearest point on the "end plane" is out of the allowed radius, find distance to the nearest edge
      	if(rval1 >0){
           point pt = pe1+Jperp2*R1;
           mindist1 = pend1.eDist(pt);
	      }else{
           point pt = pe1-Jperp2*R1;
           mindist1 = pend1.eDist(pt);
	      }
      }else{
	        //here the nearest poin in the end plane is within its radius
	        point pt = pe1+Jperp2*rval1;
	        mindist1 = pend1.eDist(pt);
      }
    }
    // for the second possibility
    if(tf2<tm1 && tf2>-tm1){
      // nearest point is potentially inside the tube
      mindist2 = jp2.eDist(pend2)-R1;
    }else{
      // project onto then "end plane" of T2
      double rval2 = (pend2-pe1).dotprod(Jperp2);
      if(std::abs(rval2) > R1){
	      // here the nearest point on the "end plane" is out of the allowed radius, find distance to the nearest edge
      	if(rval2 >0){
           point pt = pe1+Jperp2*R2;
           mindist2 = pend2.eDist(pt);
	      }else{
           point pt = pe1-Jperp2*R2;
           mindist2 = pend2.eDist(pt);
	      }
      }else{
	        //here the nearest poin in the end plane is within its radius
	        point pt = pe1+Jperp2*rval2;
	        mindist2 = pend2.eDist(pt);
      }
    }
    // update min distance (if necessary) and cheeck of we have found overlap
    if(mindist1>mindist2){
      if(mindist2<minimumDistances[inch].first[jnch-inch-1].first){
	minimumDistances[inch].first[jnch-inch-1].first =mindist2;  
      }
      if(mindist2<0.0){
        return true;
      }else{
        return false;
      }
    }else{
      if(mindist1<minimumDistances[inch].first[jnch-inch-1].first){
	minimumDistances[inch].first[jnch-inch-1].first =mindist1;  
      }
      if(mindist1<0.0){
        return true;
      }else{
        return false;
      }
    }
}



void hydrationShellMinimal::updateTruthTable(int &ind1,int &ind2,double &R){
  point T = direcList[ind1];
  point cp = midPointList[ind1];
  double tm = halfLengthList[ind1];
  for(int i=0;i<allTruthTables[ind2].size();i++){
    for(int j=0;j<allTruthTables[ind2][0].size();j++){
      point solvent = allSegments[ind2][i][j];
      point dif = solvent-cp;
      double tsol = T.dotprod(dif);
      point p = cp + T*tsol;
      double dist = p.eDist(solvent);
      if(tsol<tm && tsol>-tm && dist<= R){
	//here we have found an overlap, we add the index pair to our pair list
	if(indexList[ind2][i][j][0]==1000000){
	  indexList[ind2][i][j][0] = ind1;
	}else{
	  //check if the value is on the list
	  std::vector<int> vec =indexList[ind2][i][j];
	  if(std::find(vec.begin(), vec.end(), ind1)== vec.end()){
	  indexList[ind2][i][j].push_back(ind1);
	  }
	}
      }
      if(indexList[ind2][i][j][0]==1000000){
	allTruthTables[ind2][i][j] = 1;
      }else{
	allTruthTables[ind2][i][j] = 0;
      }
    }
  }
}



// void hydrationShellMinimal::overlapPodstruct(int i,int j,int inch,int jnch){
//   //first check if any overlap could be possible (i.e) is the centre point distnace less than the sum of the haldd length
//   point cp1 = midPointList[i];point cp2 = midPointList[j];
//   // calculate the relevant minimum distances
//   std::pair<std::vector<std::vector<double> >,std::pair<std::string,std::string> > nearAppParams;
//   point T1 = direcList[i];point T2 = direcList[j];
//   double tm1 = halfLengthList[i];double tm2 = halfLengthList[j];
//   nearAppParams =getMinTvalsFiniteNew(cp1,cp2,T1,T2,tm1,tm2);
//   bool isOverlap;
//   if(nearAppParams.second.first == "edge"){
//     bool isIn1 = checkIsIn2(cp1,cp2,T1,T2,tm1,tm2,nearAppParams.first[0][0],nearAppParams.first[0][1],inch,jnch,Rout,Rout);
//     bool isIn2 =  checkIsIn2(cp1,cp2,T1,T2,tm1,tm2,nearAppParams.first[1][0],nearAppParams.first[1][1],inch,jnch,Rout,Rout);
//     bool isIn3 = checkIsIn1(cp1,cp2,T1,T2,tm1,tm2,nearAppParams.first[2][0],nearAppParams.first[2][1],inch,jnch,Rout,Rout);
//     bool isIn4 = checkIsIn2(cp1,cp2,T1,T2,tm1,tm2,nearAppParams.first[3][0],nearAppParams.first[3][1],inch,jnch,Rout,Rout);
//     if(isIn1 == true || isIn2 ==true|| isIn3 ==true|| isIn4 ==true){
//       isOverlap = true;
//     }else{
//       isOverlap = false;
//     }
//   }else{
//     // middle overlap check
//     point p1 = cp1 + T1*nearAppParams.first[0][0];
//     point p2 = cp2 + T2*nearAppParams.first[0][1];
//     double dist = p1.eDist(p2);
//  	  minimumDistances[inch].first[jnch-inch-1].first = dist-2*Rout;  
//     if(p1.eDist(p2)< 2*Rout){
//       isOverlap = true;
//     }else{
//       isOverlap = false;
//     }
//   }
//   if(isOverlap ==true){
//     updateTruthTable(i,j,Rout);
//     updateTruthTable(j,i,Rin);
//   }
// }

void hydrationShellMinimal::overlapPodstruct(int i,int j,int inch,int jnch){
   updateTruthTable(i,j,Rout);
   updateTruthTable(j,i,Rin);
 }

void hydrationShellMinimal::overlapPodstructJoined(int i,int j){
  updateTruthTable(i,j,Rout);
  updateTruthTable(j,i,Rout);
}

void hydrationShellMinimal::allOverlap(){
  for(int i=0;i<molsize-1;i++){
    for(int j=i+1;j<molsize;j++){
      if(j==i+1){
	overlapPodstructJoined(i,j);
	//minimumDistances[i].first[0].first=-100000.0;
      }else{
	overlapPodstruct(i,j,i,j);
      }
    }
  }
  /*for(int i=0;i<molsize-1;i++){
    std::sort(minimumDistances[i].first.begin(),minimumDistances[i].first.end(),[](const std::pair<double,int> &left, const std::pair<double,int> &right) {
    return left.first < right.first;
});   
    // now find the
    //std::cout<<" sorted ? "<<i<<" "<<molsize-1<<" "<<minimumDistances[i].first.size()<<"\n";
    int k = 0;
    bool outOfRange =false;
    while(k<minimumDistances[i].first.size() && outOfRange==false){
      if(minimumDistances[i].first[k].first>0.0){
	outOfRange=true;
      }else{
	k++;
      }
    }
    minimumDistances[i].second=k;
    }*/
}

void hydrationShellMinimal::allOverlap(int &changeIndex){
  for(int i=0;i<changeIndex;i++){
    for(int j=changeIndex+1;j<molsize;j++){
      if(j==i+1){
	      overlapPodstructJoined(i,j);
      }else{
	      overlapPodstruct(i,j,i,j);  
      }
    }
  }
}


std::vector<std::pair<int,double> > hydrationShellMinimal::countDistances(int index){
  // add the zero distance components
  int zerosum = 0;
  int nps=allTruthTables[index].size();
  std::vector<std::pair<int,double> > output;
  for(int i=0;i<nps;i++){
    for(int j=0;j<ntrivs;j++){
      zerosum = zerosum + allTruthTables[index][i][j];
    }
  }
  std::pair<int,double> zp;
  zp.first = zerosum;
  zp.second = 0.0;
  output.push_back(zp);
  //now all other cases in the same layer
  if(ntrivs%2==0){
    // the even ntirvials case
    for(int l=1;l<=ntrivs/2;l++){
      int dsum=0;;
      if(l<ntrivs/2){
	for(int i=0;i<nps;i++){
	  for(int j=0;j<ntrivs;j++){
	    dsum = dsum + 2*allTruthTables[index][i][j]*allTruthTables[index][i][(j+l)%ntrivs];
	  }
	}
      }
      else{
	for(int i=0;i<nps;i++){
	  for(int j=0;j<ntrivs;j++){
	    dsum = dsum + allTruthTables[index][i][j]*allTruthTables[index][i][(j+l)%ntrivs];
	  }
	}
      }
      double dst = 2.0*Rhy*std::sin(l*Pi/ntrivs);
      std::pair<int,double> zp;
      zp.first = dsum;
      zp.second = dst;
      output.push_back(zp);
    }
  }else{
     //the odd ntirvials case
     for(int l=1;l<=(ntrivs-1)/2;l++){
       int dsum=0;
       for(int i=0;i<nps;i++){
	 for(int j=0;j<ntrivs;j++){
	   dsum = dsum + 2*allTruthTables[index][i][j]*allTruthTables[index][i][(j+l)%ntrivs];
	 }
       }
       double dst = 2.0*Rhy*std::sin(l*Pi/ntrivs);
       std::pair<int,double> zp;
       zp.first = dsum;
       zp.second = dst;
       output.push_back(zp);
    }
  }
  // Now consider the contibutions between segments only relevent if there is more than one segment;
  if(nps>1){
    // first get the vertical separation of discs;
    double cstep =cstepList[index];
    for(int k=1;k<=nps-1;k++){
      // first the cases which are directly above
      int dsum=0;
      for(int i=0;i<nps-k;i++){
	for(int j=0;j<ntrivs;j++){
	  dsum = dsum + 2*allTruthTables[index][i][j]*allTruthTables[index][i+k][j];
	}
      }
      double dst = k*cstep;
      std::pair<int,double> zp;
      zp.first = dsum;
      zp.second = dst;
      output.push_back(zp);
      //now all other cases
      if(ntrivs%2==0){
	// even ntrivial
	for(int l=1;l<=(ntrivs)/2;l++){
	  int dsum=0;
	  if(l<ntrivs/2){
	    for(int i=0;i<nps-k;i++){
	      for(int j=0;j<ntrivs;j++){
		dsum = dsum + 2*allTruthTables[index][i][j]*allTruthTables[index][i+k][(j+l)%ntrivs];
		int minin;
		if(j-l <0){
		  minin = j-l +ntrivs; 
		}else{
		  minin = j-l;
		}
		dsum = dsum + 2*allTruthTables[index][i][j]*allTruthTables[index][i+k][minin];
	      }
	    }
	  }
	  else{
	    for(int i=0;i<nps-k;i++){
	      for(int j=0;j<ntrivs;j++){
		dsum = dsum + allTruthTables[index][i][j]*allTruthTables[index][i+k][(j+l)%ntrivs];
		int minin;
		if(j-l <0){
		  minin = j -l  +ntrivs;
		}else{
		  minin = j-l;
		}
		dsum = dsum + allTruthTables[index][i][j]*allTruthTables[index][i+k][minin];
	      }
	    }
	  }
	  double dst = std::sqrt(4.0*Rhy*Rhy*std::pow(std::sin(l*Pi/ntrivs),2) + k*k*cstep*cstep);
	  std::pair<int,double> zp;
	  zp.first = dsum;
	  zp.second = dst;
	output.push_back(zp);
	}
      }else{
	for(int l=1;l<=(ntrivs-1)/2;l++){
	  int dsum=0;
	  for(int i=0;i<nps-k;i++){
	    for(int j=0;j<ntrivs;j++){
	      dsum = dsum + 2*allTruthTables[index][i][j]*allTruthTables[index][i+k][(j+l)%ntrivs];
	      int minin;
	      if(j-l <0){
		minin = j -l +ntrivs;
	      }else{
		minin = j-l;
	      }
	      dsum = dsum + 2*allTruthTables[index][i][j]*allTruthTables[index][i+k][minin];
	    }
	  }
	  double dst = std::sqrt(4.0*Rhy*Rhy*std::pow(std::sin(l*Pi/ntrivs),2) + k*k*cstep*cstep);
	  std::pair<int,double> zp;
	  zp.first = dsum;
	  zp.second = dst;
	  output.push_back(zp);
	}
      }
    }
  }
  return output;
}

double hydrationShellMinimal::helixDistFunc(int &d, int &index,int &npts){
  if(npts>1){
    double kap = mol.getCurvatureJoined(index);
    double tau = mol.getTorsionJoined(index);
    double L = mol.getAlbeadJoined(index)*npts;
    double kpt = kap*kap +tau*tau;
    double p1 = L*L*d*d*tau*tau*kpt/((npts-1)*(npts-1));
    double p2 = 2.0*kap*kap*(-1.0 + std::cos(L*d*std::sqrt(kpt)/(npts-1)));
    double denom = kpt*kpt;
    return std::sqrt((p1-p2)/denom);
  }else{
    return 0.0;
  }
}


void hydrationShellMinimal::selfScatteringSingleFast(int &nscat,double &kmin,double &kmax,int index){
  //calculate the solvent solvent distances
  double kstep = (kmax-kmin)/double(nscat-1);
  std::vector<std::pair<int,double> > solventDists = countDistances(index);
  double k;
  std::vector<std::pair<double,double> > solScatter;
  std::vector<std::pair<double,double> > molScatter;
  std::vector<std::pair<double,double> > solmolScatter;
  for(int i=0;i<nscat;i++){
    k = kmin + kstep*i;
    double sctval=0.0;
    if(k>0){
      // add distance zero values (phase is 1)
      sctval = sctval+ solventDists[0].first;
      for(int j=1;j<solventDists.size();j++){
	if(solventDists[j].first>0){
	  sctval = sctval+ solventDists[j].first*std::sin(k*solventDists[j].second)/(k*solventDists[j].second);
	}
      }
    }else{
      // k=0 all phase values are 1
      for(int j=0;j<solventDists.size();j++){
	sctval = sctval+ solventDists[j].first;
      }
    }
    std::pair<double,double> p;
    p.first = k;
    p.second = sctval;
    solScatter.push_back(p);
  }
  selfScatteringListSol[index] = solScatter;
  std::vector<point> hpts = helixpts[index];
  // now helix helix interactions; yall
  int npts =hpts.size();
  for(int i=0;i<nscat;i++){
    k = kmin + kstep*i;
    double sctval=0.0;
    if(k>0){
      sctval = npts;
      if(npts>1){
	for(int j=1;j<npts;j++){
	  double r = helixDistFunc(j,index,npts);
	  sctval = sctval + 2.0*double(npts-j)*std::sin(k*r)/(k*r);
	}
      }
    }else{
      if(npts>1){
	sctval = double(npts*npts);
      }
      else{
	sctval = npts;
      }
    }
    std::pair<double,double> p;
    p.first =k;
    p.second = sctval;
    molScatter.push_back(p);
  }
  selfScatteringListMol[index] = molScatter;
  // finally the helix solvent contirbution, I have no current spped up for this
  std::vector<point> spts;
  //make the solvent list keeping only valid solvents
  for(int i=0;i<allSegments[index].size();i++){
    for(int j=0;j<allSegments[index][0].size();j++){
      if(allTruthTables[index][i][j]==1){
	spts.push_back(allSegments[index][i][j]);
      }
    }
  }
  std::vector<double> helsoldists;
  for(int j=0;j<spts.size();j++){
    for(int i=0;i<hpts.size();i++){
      helsoldists.push_back(hpts[i].eDist(spts[j]));
    }
  }
  for(int i=0;i<nscat;i++){
    k = kmin + kstep*i;
    double sctval=0.0;
    if(k>0){
      for(int j=0;j<helsoldists.size();j++){
	if(helsoldists[j]>0){
	  sctval = sctval + 2.0*std::sin(k*helsoldists[j])/(k*helsoldists[j]);
	}else{
	  sctval = sctval + 2.0;
	}
      }
    }else{
      sctval = 2*helsoldists.size();
    }
    std::pair<double,double> p;
    p.first =k;
    p.second = sctval;
    solmolScatter.push_back(p);
  }
  selfScatteringListSolMol[index] = solmolScatter;
}

void hydrationShellMinimal::selfScatteringSingleFastRT(int &nscat,double &kmin,double &kmax,int index){
  //calculate the solvent solvent distances
  double kstep = (kmax-kmin)/double(nscat);
  std::vector<std::pair<int,double> > solventDists = countDistances(index);
  double k;
  std::vector<std::pair<double,double> > solScatter;
  std::vector<std::pair<double,double> > molScatter;
  std::vector<std::pair<double,double> > solmolScatter;
  for(int i=0;i<nscat;i++){
    k = 0.5*(kmin + kstep*i + kmin + kstep*(i+1));
    double sctval=0.0;
    if(k>0){
      // add distance zero values (phase is 1)
      sctval = sctval+ solventDists[0].first;
      for(int j=1;j<solventDists.size();j++){
      //std::cout<<"test "<<solventDists[j].first<<" "<<solventDists[j].second<<"\n";
	if(solventDists[j].first>0){
	  sctval = sctval+ solventDists[j].first*std::sin(k*solventDists[j].second)/(k*solventDists[j].second);
	}
      }
    }else{
      // k=0 all phase values are 1
      for(int j=0;j<solventDists.size();j++){
	sctval = sctval+ solventDists[j].first;
      }
    }
    std::pair<double,double> p;
    p.first = k;
    p.second = sctval;
    solScatter.push_back(p);
  }
  selfScatteringListSol[index] = solScatter;
  std::vector<point> hpts = helixpts[index];
  // now helix helix interactions; yall
  std::vector<double> heldists;
  for(int j=0;j<hpts.size()-1;j++){
    for(int i=j+1;i<hpts.size();i++){
      heldists.push_back(hpts[i].eDist(hpts[j]));
    }
  }
  for(int i=0;i<nscat;i++){
    k = 0.5*(kmin + kstep*i + kmin + kstep*(i+1));
    double sctval=0.0;
    if(k>0){
      for(int j=0;j<heldists.size();j++){
	if(heldists[j]>0.00000001){
	  sctval = sctval + 2.0*std::sin(k*heldists[j])/(k*heldists[j]);
	}else{
	  sctval = sctval + 2.0;
	}
      }
    }else{
      sctval = 2*heldists.size();
    }
    // now add the self interactions
    sctval=sctval + 2.0*double(hpts.size()); 
    std::pair<double,double> p;
    p.first =k;
    p.second = sctval;
    molScatter.push_back(p);
  }
  selfScatteringListMol[index] = molScatter;
  // finally the helix solvent contirbution, I have no current spped up for this
  std::vector<point> spts;
  //make the solvent list keeping only valid solvents
  for(int i=0;i<allSegments[index].size();i++){
    for(int j=0;j<allSegments[index][0].size();j++){
      if(allTruthTables[index][i][j]==1){
	spts.push_back(allSegments[index][i][j]);
      }
    }
  }
  std::vector<double> helsoldists;
  for(int j=0;j<spts.size();j++){
    for(int i=0;i<hpts.size();i++){
      helsoldists.push_back(hpts[i].eDist(spts[j]));
    }
  }
  for(int i=0;i<nscat;i++){
    k = 0.5*(kmin + kstep*i + kmin + kstep*(i+1));
    double sctval=0.0;
    if(k>0){
      for(int j=0;j<helsoldists.size();j++){
	if(helsoldists[j]>0){
	  sctval = sctval + 2.0*std::sin(k*helsoldists[j])/(k*helsoldists[j]);
	}else{
	  sctval = sctval + 2.0;
	}
      }
    }else{
      sctval = 2*helsoldists.size();
    }
    std::pair<double,double> p;
    p.first =k;
    p.second = sctval;
    solmolScatter.push_back(p);
  }
  selfScatteringListSolMol[index] = solmolScatter;
}

void hydrationShellMinimal::sectionAndSectionFull(int &nscat,double &kmin,double &kmax,int index1,int index2){
  int nps1 =  allSegments[index1].size();
  std::vector<point> spts1;
  //make the solvent list for section1  keeping only valid solvents
  for(int i=0;i<allSegments[index1].size();i++){
    for(int j=0;j<allSegments[index1][0].size();j++){
      if(allTruthTables[index1][i][j]==1){
	    spts1.push_back(allSegments[index1][i][j]);
      }
    }
  }
  int nps2 =  allSegments[index2].size();
  std::vector<point> spts2;
  //make the solvent list for section 2 keeping only valid solvents
  for(int i=0;i<allSegments[index2].size();i++){
    for(int j=0;j<allSegments[index2][0].size();j++){
      if(allTruthTables[index2][i][j]==1){
	    spts2.push_back(allSegments[index2][i][j]);
      }
    }
  }
  // get the helices
  std::vector<point> hpts1 = helixpts[index1];
  std::vector<point> hpts2 = helixpts[index2];
  // get the distances for helix 1 and solvents 2 
  std::vector<double> helsoldists;
  for(int j=0;j<spts2.size();j++){
    for(int i=0;i<hpts1.size();i++){
      helsoldists.push_back(hpts1[i].eDist(spts2[j]));
    }
  }
  // get the distances for helix 2 and solvents 1 
  for(int j=0;j<spts1.size();j++){
    for(int i=0;i<hpts2.size();i++){
      helsoldists.push_back(hpts2[i].eDist(spts1[j]));
    }
  }
  std::vector<double> hel1hel2dists;
  for(int j=0;j<hpts1.size();j++){
    for(int i=0;i<hpts2.size();i++){
      hel1hel2dists.push_back(hpts1[j].eDist(hpts2[i]));
    }
  }
  std::vector<double> sol1sol2dists;
  for(int j=0;j<spts1.size();j++){
    for(int i=0;i<spts2.size();i++){
      sol1sol2dists.push_back(spts1[j].eDist(spts2[i]));
    }
  }
  //helix helix scattering
  std::vector<std::pair<double,double> > helhelscat;
  double kstep = (kmax-kmin)/double(nscat-1);
  for(int i=0;i<nscat;i++){
    double k = kmin + kstep*i;
    double sctval=0.0;
    if(k>0){
      for(int j=0;j<hel1hel2dists.size();j++){
        if(hel1hel2dists[j]>0){
	         sctval = sctval + 2.0*std::sin(k*hel1hel2dists[j])/(k*hel1hel2dists[j]);
        }else{
	         sctval = sctval + 2.0;
        }
      }
    }else{
      sctval = 2*hel1hel2dists.size();
    }
    std::pair<double,double> p;
    p.first = k;
    p.second = sctval;
    helhelscat.push_back(p);
  }
  mutualScatteringListMol[index1][index2] = helhelscat;
  //solvent solvent scattering
  std::vector<std::pair<double,double> > solsolscat;
  for(int i=0;i<nscat;i++){
    double k = kmin + kstep*i;
    double sctval=0.0;
    if(k>0){
      for(int j=0;j<sol1sol2dists.size();j++){
        if(sol1sol2dists[j]>0){
	         sctval = sctval + 2.0*std::sin(k*sol1sol2dists[j])/(k*sol1sol2dists[j]);
        }else{
	         sctval = sctval + 2.0;
        }
      }
    }else{
      sctval = 2*sol1sol2dists.size();
    }
    std::pair<double,double> p;
    p.first = k;
    p.second = sctval;
    solsolscat.push_back(p);
  }
  mutualScatteringListSol[index1][index2] =solsolscat;
  std::vector<std::pair<double,double> > helsolscat;
  for(int i=0;i<nscat;i++){
    double k = kmin + kstep*i;
    double sctval=0.0;
    if(k>0){
      for(int j=0;j<helsoldists.size();j++){
        if(helsoldists[j]>0.0){
	         sctval = sctval + 2.0*std::sin(k*helsoldists[j])/(k*helsoldists[j]);
        }else{
	         sctval = sctval + 2.0;
        }
      }
    }else{
      sctval = 2*helsoldists.size();
    }
    std::pair<double,double> p;
    p.first = k;
    p.second = sctval;
    helsolscat.push_back(p);
  }
  mutualScatteringListSolMol[index1][index2] =helsolscat;
}

void hydrationShellMinimal::sectionAndSectionFullRT(int &nscat,double &kmin,double &kmax,int index1,int index2){
  int nps1 =  allSegments[index1].size();
  std::vector<point> spts1;
  //make the solvent list for section1  keeping only valid solvents
  for(int i=0;i<allSegments[index1].size();i++){
    for(int j=0;j<allSegments[index1][0].size();j++){
      if(allTruthTables[index1][i][j]==1){
	    spts1.push_back(allSegments[index1][i][j]);
      }
    }
  }
  int nps2 =  allSegments[index2].size();
  std::vector<point> spts2;
  //make the solvent list for section 2 keeping only valid solvents
  for(int i=0;i<allSegments[index2].size();i++){
    for(int j=0;j<allSegments[index2][0].size();j++){
      if(allTruthTables[index2][i][j]==1){
	    spts2.push_back(allSegments[index2][i][j]);
      }
    }
  }
  // get the helices
  std::vector<point> hpts1 = helixpts[index1];
  std::vector<point> hpts2 = helixpts[index2];
  // get the distances for helix 1 and solvents 2 
  std::vector<double> helsoldists;
  for(int j=0;j<spts2.size();j++){
    for(int i=0;i<hpts1.size();i++){
      helsoldists.push_back(hpts1[i].eDist(spts2[j]));
    }
  }
  // get the distances for helix 2 and solvents 1 
  for(int j=0;j<spts1.size();j++){
    for(int i=0;i<hpts2.size();i++){
      helsoldists.push_back(hpts2[i].eDist(spts1[j]));
    }
  }
  std::vector<double> hel1hel2dists;
  for(int j=0;j<hpts1.size();j++){
    for(int i=0;i<hpts2.size();i++){
      hel1hel2dists.push_back(hpts1[j].eDist(hpts2[i]));
    }
  }
  std::vector<double> sol1sol2dists;
  for(int j=0;j<spts1.size();j++){
    for(int i=0;i<spts2.size();i++){
      sol1sol2dists.push_back(spts1[j].eDist(spts2[i]));
    }
  }
  //helix helix scattering
  std::vector<std::pair<double,double> > helhelscat;
  double kstep = (kmax-kmin)/double(nscat);
  double k;
  for(int i=0;i<nscat;i++){
    k = 0.5*(kmin + kstep*i + kmin + kstep*(i+1));
    double sctval=0.0;
    if(k>0){
      for(int j=0;j<hel1hel2dists.size();j++){
        if(hel1hel2dists[j]>0){
	         sctval = sctval + 2.0*std::sin(k*hel1hel2dists[j])/(k*hel1hel2dists[j]);
        }else{
	         sctval = sctval + 2.0;
        }
      }
    }else{
      sctval = 2*hel1hel2dists.size();
    }
    std::pair<double,double> p;
    p.first = k;
    p.second = sctval;
    helhelscat.push_back(p);
  }
  mutualScatteringListMol[index1][index2] = helhelscat;
  //solvent solvent scattering
  std::vector<std::pair<double,double> > solsolscat;
  for(int i=0;i<nscat;i++){
    k = 0.5*(kmin + kstep*i + kmin + kstep*(i+1));
    double sctval=0.0;
    if(k>0){
      for(int j=0;j<sol1sol2dists.size();j++){
        if(sol1sol2dists[j]>0){
	         sctval = sctval + 2.0*std::sin(k*sol1sol2dists[j])/(k*sol1sol2dists[j]);
        }else{
	         sctval = sctval + 2.0;
        }
      }
    }else{
      sctval = 2*sol1sol2dists.size();
    }
    std::pair<double,double> p;
    p.first = k;
    p.second = sctval;
    solsolscat.push_back(p);
  }
  mutualScatteringListSol[index1][index2] =solsolscat;
  std::vector<std::pair<double,double> > helsolscat;
  for(int i=0;i<nscat;i++){
    k = 0.5*(kmin + kstep*i + kmin + kstep*(i+1));
    double sctval=0.0;
    if(k>0){
      for(int j=0;j<helsoldists.size();j++){
        if(helsoldists[j]>0.0){
	         sctval = sctval + 2.0*std::sin(k*helsoldists[j])/(k*helsoldists[j]);
        }else{
	         sctval = sctval + 2.0;
        }
      }
    }else{
      sctval = 2*helsoldists.size();
    }
    std::pair<double,double> p;
    p.first = k;
    p.second = sctval;
    helsolscat.push_back(p);
  }
  mutualScatteringListSolMol[index1][index2] =helsolscat;
}


void hydrationShellMinimal::sectionAndSectionUpdateLow(int &nscat,double &kmin,double &kmax,int index1,int index2){
  int nps1 =  allSegments[index1].size();
  std::vector<point> spts1;
  //make the solvent list for section1  keeping only valid solvents
  for(int i=0;i<allSegments[index1].size();i++){
    for(int j=0;j<allSegments[index1][0].size();j++){
      if(allTruthTables[index1][i][j]==1){
	    spts1.push_back(allSegments[index1][i][j]);
      }
    }
  }
  int nps2 =  allSegments[index2].size();
  std::vector<point> spts2;
  //make the solvent list for section 2 keeping only valid solvents
  for(int i=0;i<allSegments[index2].size();i++){
    for(int j=0;j<allSegments[index2][0].size();j++){
      if(allTruthTables[index2][i][j]==1){
	    spts2.push_back(allSegments[index2][i][j]);
      }
    }
  }
  // get the helices
  std::vector<point> hpts1 = helixpts[index1];
  std::vector<point> hpts2 = helixpts[index2];
  // get the distances for helix 1 and solvents 2 
  std::vector<double> helsoldists;
  for(int j=0;j<spts2.size();j++){
    for(int i=0;i<hpts1.size();i++){
      helsoldists.push_back(hpts1[i].eDist(spts2[j]));
    }
  }
  // get the distances for helix 2 and solvents 1 
  for(int j=0;j<spts1.size();j++){
    for(int i=0;i<hpts2.size();i++){
      helsoldists.push_back(hpts2[i].eDist(spts1[j]));
    }
  }
  std::vector<double> hel1hel2dists;
  for(int j=0;j<hpts1.size();j++){
    for(int i=0;i<hpts2.size();i++){
      hel1hel2dists.push_back(hpts1[j].eDist(hpts2[i]));
    }
  }
  std::vector<double> sol1sol2dists;
  for(int j=0;j<spts1.size();j++){
    for(int i=0;i<spts2.size();i++){
      sol1sol2dists.push_back(spts1[j].eDist(spts2[i]));
    }
  }
  double kstep = (kmax-kmin)/(nscat-1);
  //solvent solvent scattering
  std::vector<std::pair<double,double> > solsolscat;
  for(int i=0;i<nscat;i++){
    double k = kmin + kstep*i;
    double sctval=0.0;
    if(k>0){
      for(int j=0;j<sol1sol2dists.size();j++){
        if(sol1sol2dists[j]>0){
	         sctval = sctval + 2.0*std::sin(k*sol1sol2dists[j])/(k*sol1sol2dists[j]);
        }else{
	         sctval = sctval + 2.0;
        }
      }
    }else{
      sctval = 2*sol1sol2dists.size();
    }
    std::pair<double,double> p;
    p.first = k;
    p.second = sctval;
    solsolscat.push_back(p);
  }
  mutualScatteringListSol[index1][index2] =solsolscat;
  std::vector<std::pair<double,double> > helsolscat;
  for(int i=0;i<nscat;i++){
    double k = kmin + kstep*i;
    double sctval=0.0;
    if(k>0){
      for(int j=0;j<helsoldists.size();j++){
        if(helsoldists[j]>0){
	         sctval = sctval + 2.0*std::sin(k*helsoldists[j])/(k*helsoldists[j]);
        }else{
	         sctval = sctval + 2.0;
        }
      }
    }else{
      sctval = 2*helsoldists.size();
    }
    std::pair<double,double> p;
    p.first = k;
    p.second = sctval;
    helsolscat.push_back(p);
  }
  mutualScatteringListSolMol[index1][index2] =helsolscat;
}

void hydrationShellMinimal::sectionAndSectionUpdateLowRT(int &nscat,double &kmin,double &kmax,int index1,int index2){
  int nps1 =  allSegments[index1].size();
  std::vector<point> spts1;
  //make the solvent list for section1  keeping only valid solvents
  for(int i=0;i<allSegments[index1].size();i++){
    for(int j=0;j<allSegments[index1][0].size();j++){
      if(allTruthTables[index1][i][j]==1){
	    spts1.push_back(allSegments[index1][i][j]);
      }
    }
  }
  int nps2 =  allSegments[index2].size();
  std::vector<point> spts2;
  //make the solvent list for section 2 keeping only valid solvents
  for(int i=0;i<allSegments[index2].size();i++){
    for(int j=0;j<allSegments[index2][0].size();j++){
      if(allTruthTables[index2][i][j]==1){
	    spts2.push_back(allSegments[index2][i][j]);
      }
    }
  }
  // get the helices
  std::vector<point> hpts1 = helixpts[index1];
  std::vector<point> hpts2 = helixpts[index2];
  // get the distances for helix 1 and solvents 2 
  std::vector<double> helsoldists;
  for(int j=0;j<spts2.size();j++){
    for(int i=0;i<hpts1.size();i++){
      helsoldists.push_back(hpts1[i].eDist(spts2[j]));
    }
  }
  // get the distances for helix 2 and solvents 1 
  for(int j=0;j<spts1.size();j++){
    for(int i=0;i<hpts2.size();i++){
      helsoldists.push_back(hpts2[i].eDist(spts1[j]));
    }
  }
  std::vector<double> hel1hel2dists;
  for(int j=0;j<hpts1.size();j++){
    for(int i=0;i<hpts2.size();i++){
      hel1hel2dists.push_back(hpts1[j].eDist(hpts2[i]));
    }
  }
  std::vector<double> sol1sol2dists;
  for(int j=0;j<spts1.size();j++){
    for(int i=0;i<spts2.size();i++){
      sol1sol2dists.push_back(spts1[j].eDist(spts2[i]));
    }
  }
  double kstep = (kmax-kmin)/(nscat);
  //solvent solvent scattering
  std::vector<std::pair<double,double> > solsolscat;
  double k;
  for(int i=0;i<nscat;i++){
    k = 0.5*(kmin + kstep*i + kmin + kstep*(i+1));
    double sctval=0.0;
    if(k>0){
      for(int j=0;j<sol1sol2dists.size();j++){
        if(sol1sol2dists[j]>0){
	         sctval = sctval + 2.0*std::sin(k*sol1sol2dists[j])/(k*sol1sol2dists[j]);
        }else{
	         sctval = sctval + 2.0;
        }
      }
    }else{
      sctval = 2*sol1sol2dists.size();
    }
    std::pair<double,double> p;
    p.first = k;
    p.second = sctval;
    solsolscat.push_back(p);
  }
  mutualScatteringListSol[index1][index2] =solsolscat;
  std::vector<std::pair<double,double> > helsolscat;
  for(int i=0;i<nscat;i++){
    k = 0.5*(kmin + kstep*i + kmin + kstep*(i+1));
    double sctval=0.0;
    if(k>0){
      for(int j=0;j<helsoldists.size();j++){
        if(helsoldists[j]>0){
	         sctval = sctval + 2.0*std::sin(k*helsoldists[j])/(k*helsoldists[j]);
        }else{
	         sctval = sctval + 2.0;
        }
      }
    }else{
      sctval = 2*helsoldists.size();
    }
    std::pair<double,double> p;
    p.first = k;
    p.second = sctval;
    helsolscat.push_back(p);
  }
  mutualScatteringListSolMol[index1][index2] =helsolscat;
}



void hydrationShellMinimal::selfScatteringAll(int &nscat,double &kmin,double &kmax){
  for(int i=0;i<molsize;i++){
    selfScatteringSingleFast(nscat,kmin,kmax,i);
    addOneScatterSelf(nscat,i);
  }
}

void hydrationShellMinimal::mutualScatteringAll(int &nscat,double &kmin,double &kmax,int index){
  for(int i = index+1;i<molsize;i++){
    sectionAndSectionFull(nscat,kmin,kmax,index,i);
    addOneScatterMutual(nscat,index,i);
  }
}

void hydrationShellMinimal::selfScatteringAllRT(int &nscat,double &kmin,double &kmax){
  for(int i=0;i<molsize;i++){
    selfScatteringSingleFastRT(nscat,kmin,kmax,i);
    addOneScatterSelf(nscat,i);
  }
}

void hydrationShellMinimal::mutualScatteringAllRT(int &nscat,double &kmin,double &kmax,int index){
  for(int i = index+1;i<molsize;i++){
    sectionAndSectionFullRT(nscat,kmin,kmax,index,i);
    addOneScatterMutual(nscat,index,i);
  }
}


void hydrationShellMinimal::getAllScatter(int &nscat,double &kmin,double &kmax){
  getAllHelices();
  // we should ignore some of the linker scattering form the secons just before the  (specifically its helices
  double kstep = (kmax-kmin)/double(nscat-1);
  for(int i =0;i< nscat;i++){
    double k = kmin +  kstep*i;
    std::pair<double,double> p;
    p.first = k;
    p.second =0.0;
    totalScatterSol.push_back(p);
    totalScatterMol.push_back(p);
    totalScatterSolMol.push_back(p);
    emptyScatterSol.push_back(p);
    emptyScatterMol.push_back(p);
    emptyScatterSolMol.push_back(p);
  }
  //std::vector<int> ignorePts = mol.getIgnorePts();		  
  selfScatteringAll(nscat,kmin,kmax);
  //std::cout<<"done self scattering\n";
  for(int i=0;i<molsize-1;i++){
    mutualScatteringAll(nscat,kmin,kmax,i);
  }
}

void hydrationShellMinimal::getAllScatterRT(int &nscat,double &kmin,double &kmax){
  getAllHelices();
  // we should ignore some of the linker scattering form the secons just before the  (specifically its helices
  double kstep = (kmax-kmin)/double(nscat);
  for(int i =0;i< nscat;i++){
    double k = 0.5*(kmin + kstep*i + kmin + kstep*(i+1));
    std::pair<double,double> p;
    p.first = k;
    p.second =0.0;
    totalScatterSol.push_back(p);
    totalScatterMol.push_back(p);
    totalScatterSolMol.push_back(p);
    emptyScatterSol.push_back(p);
    emptyScatterMol.push_back(p);
    emptyScatterSolMol.push_back(p);
  }
  //std::vector<int> ignorePts = mol.getIgnorePts();		  
  selfScatteringAllRT(nscat,kmin,kmax);
  //std::cout<<"done self scattering\n";
  for(int i=0;i<molsize-1;i++){
    mutualScatteringAllRT(nscat,kmin,kmax,i);
  }
}



void hydrationShellMinimal::getInitialScatter(int &nk,double &kmin,double &kmax){
  tubeParamList();
  constructInitialState();
  allOverlap();
  getAllScatter(nk,kmin,kmax);
}

void hydrationShellMinimal::getInitialScatterRT(int &nk,double &kmin,double &kmax){
  tubeParamList();
  constructInitialState();
  allOverlap();
  getAllScatterRT(nk,kmin,kmax);
}


void hydrationShellMinimal::solventMoleculeDistances(std::vector<double> &molSolDistances,std::vector<double> &solSolDistances){
  // first find all solvent moelcule distances
  for(int k=0;k<allSegments.size();k++){
    for(int i=0;i<allSegments[k].size();i++){
      for(int j=0;j<allSegments[k][0].size();j++){
	bool overlapped = false;
	std::vector<double> solMolPosDists;
	for(int l=0;l<helixpts.size();l++){
	  for(int m=0;m<helixpts[l].size();m++){
	    // declare molpt
	    if(k!=l){
	      double dist =helixpts[l][m].eDist(allSegments[k][i][j]);
	      if(dist<5.4){
		overlapped=true;
		//std::cout<<dist<<" "<<k<<" "<<i<<" "<<j<<" "<<l<<" "<<m<<"\n"; 
	      }else{
		solMolPosDists.push_back(dist);
	      }
	    }else{
	       double dist =helixpts[l][m].eDist(allSegments[k][i][j]);
	       solMolPosDists.push_back(dist);
	    }
	  }
	}
	if(overlapped==true){
	  allTruthTables[k][i][j]=0;
	}else{
	  // this solvent is still in play, store its distances.
	  if(molSolDistances.size()==0){
	    molSolDistances= solMolPosDists;
	  }else{
	    molSolDistances.insert(molSolDistances.end(),solMolPosDists.begin(),solMolPosDists.end());
	  }
	}
      }
    }
  }
  for(int k=0;k<allSegments.size()-1;k++){
    for(int i=0;i<allSegments[k].size();i++){
      for(int j=0;j<allSegments[k][0].size();j++){
	//std::cout<<allTruthTables[k][i][j]<<"\n";
	if(allTruthTables[k][i][j]==1){
	  for(int l=k+1;l<allSegments.size();l++){
	    for(int m=0;m<allSegments[l].size();m++){
	      for(int n=0;n<allSegments[l][0].size();n++){
		if(allTruthTables[l][m][n]==1){
		 double soldist = allSegments[l][m][n].eDist(allSegments[k][i][j]);
		 //std::cout<<soldist<<"\n";
		 if(soldist<2.5){
		   allTruthTables[l][m][n]=0;
		 }else{
		   solSolDistances.push_back(soldist);
		 }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  //std::cout<<"escape ? "<<molSolDistances.size()<<" "<<solSolDistances.size()<<"\n";
  // now we check for solvent solvent distances *between* sections (not within a section)
}



std::vector<std::vector<point> > hydrationShellMinimal::returnFlatSolList(){
  solPtsFlat.clear();
  int noflps =0;
  for(int k=0;k<allSegments.size();k++){
    std::vector<point> spts;
    for(int i=0;i<allSegments[k].size();i++){
      for(int j=0;j<allSegments[k][0].size();j++){
	//std::cout<<i<<" "<<j<<" "<<k<<" "<<allTruthTables[k][i][j]<<"\n";
	if(allTruthTables[k][i][j]==1){
	  spts.push_back(allSegments[k][i][j]);
	}
      }
    }
    solPtsFlat.push_back(spts);
    //noflps = noflps + spts.size();
    //std::cout<<"in full segment "<<k<<" "<<spts.size()<<"\n";
      //also add the helical molecules to the flat helix list
  }
  // std::cout<<"full calc "<<noflps<<"\n";
  return solPtsFlat;
}

void hydrationShellMinimal::getFlatSolList(){
  solPtsFlat.clear();
  helPtsFlat.clear();
  int noflps =0;
  for(int k=0;k<allSegments.size();k++){
    std::vector<point> spts;
    for(int i=0;i<allSegments[k].size();i++){
      for(int j=0;j<allSegments[k][0].size();j++){
	if(allTruthTables[k][i][j]==1){
	  spts.push_back(allSegments[k][i][j]);
	}
      }
    }
    solPtsFlat.push_back(spts);
    //noflps = noflps + spts.size();
    //std::cout<<"in full segment "<<k<<" "<<spts.size()<<"\n";
      //also add the helical molecules to the flat helix list
    helPtsFlat.push_back(helixpts[k]);
  }
  // std::cout<<"full calc "<<noflps<<"\n";
}

std::vector<double> hydrationShellMinimal::flattenVector(std::vector<std::vector<std::vector<double> > > &vec){
  std::vector<double> vecOut;
  vecOut = vec[0][0];
  for(int j=1;j<vec[0].size();j++){
      vecOut.insert(vecOut.end(),vec[0][j].begin(),vec[0][j].end());
  }
  for(int i=1;i<vec.size();i++){
    for(int j=0;j<vec[i].size();j++){
      vecOut.insert(vecOut.end(),vec[i][j].begin(),vec[i][j].end());
    }
  }
  return vecOut;
}

std::vector<int> hydrationShellMinimal::binDistancesMolecule(){
  isTooClose=false;
  std::vector<double> distances;
  std::vector<double> subDistances;
  std::vector<std::vector<double> > subDistanceSet;
  molmolDistSum.clear();
  for(int i=0;i<helPtsFlat.size()-1;i++){
    for(int j=i+1;j<helPtsFlat.size();j++){
      for(int k=0;k<helPtsFlat[i].size();k++){
	      for(int l=0;l<helPtsFlat[j].size();l++){
          double dst =   helPtsFlat[i][k].eDist(helPtsFlat[j][l]);
	        subDistances.push_back(dst);
         if(j>i+1 && dst < mutualDistCutOff){
            //std::cout<<i<<" "<<j<<" "<<dst<<" "<<mutualDistCutOff<<"\n";
            isTooClose=true;
          }
        }
      }
      subDistanceSet.push_back(subDistances);
      subDistances.clear();
    }
    molmolDistSum.push_back(subDistanceSet);
    subDistanceSet.clear();
  }
  distances = flattenVector(molmolDistSum);
  std::sort(distances.begin(),distances.end());
  double binsize = 0.2;
  int nbins = int(ceil(maxdist/0.2));
  std::vector<double> binlist;
  for(int i =0;i<=nbins;i++){
    binlist.push_back(binsize*i);
  }
  int k=0;int kprev=0;
  int dlen = distances.size();
  std::vector<int> binSizes;
  for(int i=0;i<nbins;i++){
    while(distances[k]<=binlist[i+1] && k<dlen){
      k++;
    }
    if((k-kprev)>=1){
      int noInBin = k-kprev;
      binSizes.push_back(noInBin);
      molIntbound =i;
    }else{
      binSizes.push_back(0);
    }
    kprev = k;
  }
  return binSizes;
}

bool hydrationShellMinimal::isCutOffViolated(){
  return isTooClose;
}

double hydrationShellMinimal::getCuttOffDistMutual(){
  return mutualDistCutOff;
}

std::vector<int> hydrationShellMinimal::binDistancesSolvent(){
  std::vector<double> distances;
  std::vector<double> subDistances;
  std::vector<std::vector<double> > subDistanceSet;
  solsolDistSum.clear();
  for(int i=0;i<solPtsFlat.size()-1;i++){
    for(int j=i+1;j<solPtsFlat.size();j++){
      for(int k=0;k<solPtsFlat[i].size();k++){
	for(int l=0;l<solPtsFlat[j].size();l++){
	  subDistances.push_back(solPtsFlat[i][k].eDist(solPtsFlat[j][l]));
	}
      }
      subDistanceSet.push_back(subDistances);
      subDistances.clear();
    }
    solsolDistSum.push_back(subDistanceSet);
    subDistanceSet.clear();
  }
  distances = flattenVector(solsolDistSum);
  std::sort(distances.begin(),distances.end());
  double binsize = 0.2;
  int nbins = int(ceil(maxdist/0.2));
  std::vector<double> binlist;
  for(int i =0;i<=nbins;i++){
    binlist.push_back(binsize*i);
  }
  int k=0;int kprev=0;
  int dlen = distances.size();
  std::vector<int> binSizes;
  for(int i=0;i<nbins;i++){
    while(distances[k]<=binlist[i+1] && k<dlen){
      k++;
    }
    if((k-kprev)>=1){
      int noInBin = k-kprev;
      binSizes.push_back(noInBin);
      solIntbound =i;
    }else{
      binSizes.push_back(0);
    }
    kprev = k;
  }
  return binSizes;
}

std::vector<int> hydrationShellMinimal::binDistancesSolventMolecule(){
  std::vector<double> distances2;
  std::vector<double> distances;
  std::vector<double> subDistances;
  std::vector<std::vector<double> > subDistanceSet;
  molsolDistSum.clear();
  solmolDistSum.clear();
  for(int i=0;i<helPtsFlat.size()-1;i++){
    for(int j=i+1;j<solPtsFlat.size();j++){
      for(int k=0;k<helPtsFlat[i].size();k++){
	for(int l=0;l<solPtsFlat[j].size();l++){
	  subDistances.push_back(helPtsFlat[i][k].eDist(solPtsFlat[j][l]));
	}
      }
      subDistanceSet.push_back(subDistances);
      subDistances.clear();
    }
    molsolDistSum.push_back(subDistanceSet);
    subDistanceSet.clear();
  }
  for(int i=0;i<solPtsFlat.size()-1;i++){
    for(int j=i+1;j<helPtsFlat.size();j++){
      int npts =0;
      for(int k=0;k<solPtsFlat[i].size();k++){
	for(int l=0;l<helPtsFlat[j].size();l++){
	  subDistances.push_back(solPtsFlat[i][k].eDist(helPtsFlat[j][l]));
	}
      }
      subDistanceSet.push_back(subDistances);
      subDistances.clear();
    }
    solmolDistSum.push_back(subDistanceSet);
    subDistanceSet.clear();
  }
  distances = flattenVector(molsolDistSum);
  distances2 = flattenVector(solmolDistSum);
  distances.insert(distances.end(),distances2.begin(),distances2.end());
  std::sort(distances.begin(),distances.end());
  double binsize = 0.2;
  int nbins = int(ceil(maxdist/0.2));
  std::vector<double> binlist;
  for(int i =0;i<=nbins;i++){
    binlist.push_back(binsize*i);
  }
  int k=0;int kprev=0;
  int dlen = distances.size();
  std::vector<int> binSizes;
  for(int i=0;i<nbins;i++){
    while(distances[k]<=binlist[i+1] && k<dlen){
      k++;
    }
    if((k-kprev)>=1){
      int noInBin = k-kprev;
      binSizes.push_back(noInBin);
      solmolIntbound =i;
    }else{
      binSizes.push_back(0);
    }
    kprev = k;
  }
  return binSizes;
}

void hydrationShellMinimal::setNScat(double &kmin,double &kmax){
  nscat = std::round(1.0*maxScatDist*(kmax-kmin)*0.31830988618);
}

void hydrationShellMinimal::preCalculatePhases(double &kmin,double &kmax){
  double kstep = (kmax-kmin)/double(nscat);
  double binsize = 0.2;
  int nbins = int(ceil(maxdist/0.2));
  std::vector<double> binlist;
  double kprod,k,dist;
  phaseList.resize(nscat,std::vector<double>(nbins,0.0));
  for(int j=0;j<nscat;j++){
    k = kmin + kstep*j;
    if(k>0){
      for(int i=0;i<nbins;i++){
	      dist = 0.2*i+0.1;
	      kprod = k*dist;
	      phaseList[j][i] = 2.0*(std::sin(kprod)/kprod);
      }
    }else{
      for(int i=0;i<nbins;i++){
	    phaseList[j][i] = 2.0;
      }
    }
  }
}

void hydrationShellMinimal::preCalculatePhasesRT(double &kmin,double &kmax){
  double kstep = (kmax-kmin)/double(nscat);
  double binsize = 0.2;
  int nbins = int(ceil(maxdist/0.2));
  std::vector<double> binlist;
  double kprod,k,dist;
  phaseList.resize(nscat,std::vector<double>(nbins,0.0));
  for(int j=0;j<nscat;j++){
    k = 0.5*(kmin + kstep*j + kmin + kstep*(j+1));
    if(k>0){
      for(int i=0;i<nbins;i++){
	dist = 0.2*i+0.1;
	kprod = k*dist;
	phaseList[j][i] = 2.0*(std::sin(kprod)/kprod);
      }
    }else{
      for(int i=0;i<nbins;i++){
	phaseList[j][i] = 2.0;
      }
    }
  }
}


void hydrationShellMinimal::addScatteringLists(std::vector<std::pair<double,double> > &l1,std::vector<std::pair<double,double> > &l2){
  //here we add list l2 to l1
  for(int i=0;i<l1.size();i++){
    l1[i].second = l1[i].second + l2[i].second;
  }
}



void hydrationShellMinimal::getScatteringBinnedRT(double &kmin,double &kmax){
  tubeParamList();
  constructInitialState();
  allOverlap();
  getAllHelices();
  getFlatSolList();
  //now the mutual scattering, first fill the distance bins
  distBinsMol = binDistancesMolecule();
  distBinsSol = binDistancesSolvent();
  distBinsSolMol = binDistancesSolventMolecule();
  kminVal = kmin;
  kmaxVal =kmax;
  //nscat = std::round(2.0*maxdist*(kmax-kmin)*0.31830988618);
  double kstep = (kmax-kmin)/double(nscat);
  for(int i =0;i< nscat;i++){
    double k = 0.5*(kmin + kstep*i + kmin + kstep*(i+1));
    std::pair<double,double> p;
    p.first = k;
    p.second =0.0;
    totalScatterSol.push_back(p);
    totalScatterMol.push_back(p);
    totalScatterSolMol.push_back(p);
    emptyScatterSol.push_back(p);
    emptyScatterMol.push_back(p);
    emptyScatterSolMol.push_back(p);
  }
  //preCalculatePhasesRT(nscat,kmin,kmax);
  selfScatteringAllRT(nscat,kmin,kmax);
  std::pair<double,double> scatPair;
  mutualScatterMol.clear();
  mutualScatterSol.clear();
  mutualScatterSolMol.clear();
  // now calculate the scattering
  for(int j=0;j<nscat;j++){
    double scatVal =0.0;
    for(int i=0;i<=molIntbound;i++){
      scatVal = scatVal + phaseList[j][i]*distBinsMol[i];
    }
    scatPair.first = 0.5*(kmin + kstep*j + kmin + kstep*(j+1));
    scatPair.second = scatVal;
    mutualScatterMol.push_back(scatPair);
  }
  for(int j=0;j<nscat;j++){
    double scatVal =0.0;
    for(int i=0;i<=solIntbound;i++){
      scatVal = scatVal + phaseList[j][i]*distBinsSol[i];
    }
    scatPair.first = 0.5*(kmin + kstep*j + kmin + kstep*(j+1));
    scatPair.second = scatVal;
    mutualScatterSol.push_back(scatPair);
  }
  for(int j=0;j<nscat;j++){
    double scatVal =0.0;
    for(int i=0;i<solmolIntbound;i++){
      scatVal = scatVal + phaseList[j][i]*distBinsSolMol[i];
    }
    scatPair.first = 0.5*(kmin + kstep*j + kmin + kstep*(j+1));
    scatPair.second = scatVal;
    mutualScatterSolMol.push_back(scatPair);
  }
  addScatteringLists(totalScatterMol,mutualScatterMol);
  addScatteringLists(totalScatterSol,mutualScatterSol);
  addScatteringLists(totalScatterSolMol,mutualScatterSolMol);
}

void hydrationShellMinimal::getScatteringBinnedRTUpdate(){
  //set up
  tubeParamList();
  constructInitialState();
  allOverlap();
  getAllHelices();
  getFlatSolList();
  double kmin = kminVal;
  double kmax = kmaxVal;
  //now the mutual scattering, first fill the distance bins
  distBinsMol = binDistancesMolecule();
  distBinsSol = binDistancesSolvent();
  distBinsSolMol = binDistancesSolventMolecule();
  //clear previous scattering vectors
  //nscat = std::round(2.0*maxdist*(kmax-kmin)*0.31830988618);
  double kstep = (kmax-kmin)/double(nscat);
  totalScatterSol.clear();
  totalScatterMol.clear();
  totalScatterSolMol.clear();
  emptyScatterSol.clear();
  emptyScatterMol.clear();
  emptyScatterSolMol.clear();
  mutualScatterMol.clear();
  mutualScatterSol.clear();
  mutualScatterSolMol.clear();
  for(int i =0;i< nscat;i++){
    double k = 0.5*(kmin + kstep*i + kmin + kstep*(i+1));
    std::pair<double,double> p;
    p.first = k;
    p.second =0.0;
    totalScatterSol.push_back(p);
    totalScatterMol.push_back(p);
    totalScatterSolMol.push_back(p);
    emptyScatterSol.push_back(p);
    emptyScatterMol.push_back(p);
    emptyScatterSolMol.push_back(p);
  }
  //preCalculatePhasesRT(nscat,kmin,kmax);
  selfScatteringAllRT(nscat,kmin,kmax);
  std::pair<double,double> scatPair;
  // now calculate the scattering
  /*for(int i=0;i<totalScatterSol.size();i++){
    std::cout<<" old "<<totalScatterSol[i].first<<" "<<totalScatterSol[i].second<<"\n";
  }*/
  for(int j=0;j<nscat;j++){
    double scatVal =0.0;
    for(int i=0;i<=molIntbound;i++){
      scatVal = scatVal + phaseList[j][i]*distBinsMol[i];
    }
    scatPair.first = 0.5*(kmin + kstep*j + kmin + kstep*(j+1));
    scatPair.second = scatVal;
    mutualScatterMol.push_back(scatPair);
  }
  for(int j=0;j<nscat;j++){
    double scatVal =0.0;
    for(int i=0;i<=solIntbound;i++){
      scatVal = scatVal + phaseList[j][i]*distBinsSol[i];
    }
    scatPair.first = 0.5*(kmin + kstep*j + kmin + kstep*(j+1));
    scatPair.second = scatVal;
    mutualScatterSol.push_back(scatPair);
  }
  for(int j=0;j<nscat;j++){
    double scatVal =0.0;
    for(int i=0;i<solmolIntbound;i++){
      scatVal = scatVal + phaseList[j][i]*distBinsSolMol[i];
    }
    scatPair.first = 0.5*(kmin + kstep*j + kmin + kstep*(j+1));
    scatPair.second = scatVal;
    mutualScatterSolMol.push_back(scatPair);
  }
  addScatteringLists(totalScatterMol,mutualScatterMol);
  addScatteringLists(totalScatterSol,mutualScatterSol);
  addScatteringLists(totalScatterSolMol,mutualScatterSolMol);
  /*std::cout<<"sol"<<"\n";
  for(int i =0;i<totalScatterSol.size();i++){
    std::cout<<"TST "<<totalScatterSol[i].first<<" "<<totalScatterSol[i].second<<"\n";
  }
  std::cout<<"mol"<<"\n";
  for(int i =0;i<totalScatterMol.size();i++){
    std::cout<<"TST "<<totalScatterSol[i].first<<" "<<totalScatterMol[i].second<<"\n";
  }
  std::cout<<"solmol"<<"\n";
  for(int i =0;i<totalScatterSolMol.size();i++){
    std::cout<<"TST "<<totalScatterSol[i].first<<" "<<totalScatterSolMol[i].second<<"\n";
  }*/
}

void hydrationShellMinimal::getScatteringBinnedRTUpdate(int &changeIndex){
  //set up
  tubeParamList();
  constructInitialState();
  allOverlap(changeIndex);
  getAllHelices();
  getFlatSolList();
  double kmin = kminVal;
  double kmax = kmaxVal;
  //now the mutual scattering, first fill the distance bins
  distBinsMol = binDistancesMolecule();
  distBinsSol = binDistancesSolvent();
  distBinsSolMol = binDistancesSolventMolecule();
  //clear previous scattering vectors
  //nscat = std::round(2.0*maxdist*(kmax-kmin)*0.31830988618);
  double kstep = (kmax-kmin)/double(nscat);
  totalScatterSol.clear();
  totalScatterMol.clear();
  totalScatterSolMol.clear();
  emptyScatterSol.clear();
  emptyScatterMol.clear();
  emptyScatterSolMol.clear();
  mutualScatterMol.clear();
  mutualScatterSol.clear();
  mutualScatterSolMol.clear();
  for(int i =0;i< nscat;i++){
    double k = 0.5*(kmin + kstep*i + kmin + kstep*(i+1));
    std::pair<double,double> p;
    p.first = k;
    p.second =0.0;
    totalScatterSol.push_back(p);
    totalScatterMol.push_back(p);
    totalScatterSolMol.push_back(p);
    emptyScatterSol.push_back(p);
    emptyScatterMol.push_back(p);
    emptyScatterSolMol.push_back(p);
  }
  //preCalculatePhasesRT(nscat,kmin,kmax);
  selfScatteringAllRT(nscat,kmin,kmax);
  std::pair<double,double> scatPair;
  // now calculate the scattering
  /*for(int i=0;i<totalScatterSol.size();i++){
    std::cout<<" old "<<totalScatterSol[i].first<<" "<<totalScatterSol[i].second<<"\n";
  }*/
  for(int j=0;j<nscat;j++){
    double scatVal =0.0;
    for(int i=0;i<=molIntbound;i++){
      scatVal = scatVal + phaseList[j][i]*distBinsMol[i];
    }
    scatPair.first = 0.5*(kmin + kstep*j + kmin + kstep*(j+1));
    scatPair.second = scatVal;
    mutualScatterMol.push_back(scatPair);
  }
  for(int j=0;j<nscat;j++){
    double scatVal =0.0;
    for(int i=0;i<=solIntbound;i++){
      scatVal = scatVal + phaseList[j][i]*distBinsSol[i];
    }
    scatPair.first = 0.5*(kmin + kstep*j + kmin + kstep*(j+1));
    scatPair.second = scatVal;
    mutualScatterSol.push_back(scatPair);
  }
  for(int j=0;j<nscat;j++){
    double scatVal =0.0;
    for(int i=0;i<solmolIntbound;i++){
      scatVal = scatVal + phaseList[j][i]*distBinsSolMol[i];
    }
    scatPair.first = 0.5*(kmin + kstep*j + kmin + kstep*(j+1));
    scatPair.second = scatVal;
    mutualScatterSolMol.push_back(scatPair);
  }
  addScatteringLists(totalScatterMol,mutualScatterMol);
  addScatteringLists(totalScatterSol,mutualScatterSol);
  addScatteringLists(totalScatterSolMol,mutualScatterSolMol);
  /*std::cout<<"sol"<<"\n";
  for(int i =0;i<totalScatterSol.size();i++){
    std::cout<<"TST "<<totalScatterSol[i].first<<" "<<totalScatterSol[i].second<<"\n";
  }
  std::cout<<"mol"<<"\n";
  for(int i =0;i<totalScatterMol.size();i++){
    std::cout<<"TST "<<totalScatterSol[i].first<<" "<<totalScatterMol[i].second<<"\n";
  }
  std::cout<<"solmol"<<"\n";
  for(int i =0;i<totalScatterSolMol.size();i++){
    std::cout<<"TST "<<totalScatterSol[i].first<<" "<<totalScatterSolMol[i].second<<"\n";
  }*/
}




void hydrationShellMinimal::addOneScatterSelf(int &nscat,int &index){
    for(int j=0;j<nscat;j++){
      totalScatterSol[j].second =totalScatterSol[j].second+selfScatteringListSol[index][j].second;
      totalScatterMol[j].second =totalScatterMol[j].second+selfScatteringListMol[index][j].second;
      totalScatterSolMol[j].second =totalScatterSolMol[j].second+selfScatteringListSolMol[index][j].second;
  }
}

void hydrationShellMinimal::addOneScatterMutual(int &nscat,int &index1,int &index2){
   for(int k=0;k<nscat;k++){
     totalScatterSol[k].second =totalScatterSol[k].second+mutualScatteringListSol[index1][index2][k].second;
     totalScatterMol[k].second =totalScatterMol[k].second+mutualScatteringListMol[index1][index2][k].second;
     totalScatterSolMol[k].second =totalScatterSolMol[k].second+mutualScatteringListSolMol[index1][index2][k].second;
   }
}



void hydrationShellMinimal::printSelfDistances(int index){
  std::vector<std::pair<int,double> > tdists = countDistances(index);
  for(int i=0;i<tdists.size();i++){
    std::pair<int,double> tp = tdists[i];
    std::cout<<tp.first<<" "<<tp.second<<"\n";
  }
}

void hydrationShellMinimal::printOneSelfScatter(int index){
  std::vector<std::pair<double,double> > sssol = selfScatteringListMol[index];
  for(int i=0;i<sssol.size();i++){
    std::pair<double,double> p = sssol[i];
      std::cout<<p.first<<" "<<p.second<<"\n";
  }
  /*std::cout<<"self scatter section "<<index<<"solvent with Molecule:\n";
  std::vector<std::pair<double,double> > sssolmol = selfScatteringListSolMol[index];
  for(int i=0;i<sssolmol.size();i++){
    std::pair<double,double> p = sssolmol[i];
    std::cout<<p.first<<" "<<p.second<<"\n";
  }
  std::vector<std::pair<double,double> > ssmol = selfScatteringListSol[index];
  std::cout<<"self scatter section "<<index<<"molecule:\n";
  for(int i=0;i<ssmol.size();i++){
    std::pair<double,double> p = ssmol[i];
    std::cout<<p.first<<" "<<p.second<<"\n";
    }*/
}


std::vector<std::pair<double,double> > hydrationShellMinimal::getTotalSelfScatter(){
  double soltot=0.0;
  double solmoltot=0.0;
  double moltot=0.0;
  std::cout<<"in here ?\n";
  std::vector<std::pair<double,double> > outList;
  for(int i=0;i<selfScatteringListMol.size();i++){
    std::vector<std::pair<double,double> > sssol = selfScatteringListMol[i];
    for(int j=0;j<sssol.size();j++){
      if(i>0){
	std::pair<double,double> p = sssol[j];
	outList[j].second = outList[j].second + p.second;
      }else{
	outList = sssol;
      }
    }
  }
  return outList;
}


void hydrationShellMinimal::printOneMutualScatter(int index1,int index2){
  std::vector<std::pair<double,double> > sssol = mutualScatteringListMol[index1][index2];
  std::cout<<"mutual scatter section "<<index1<<" "<<index2<<"solvent:\n";
  for(int i=0;i<sssol.size();i++){
    std::pair<double,double> p = sssol[i];
    std::cout<<p.first<<" "<<p.second<<"\n";
  }
  std::cout<<"mutual scatter section "<<index1<<" "<<index2<<"solvent with Molecule:\n";
  std::vector<std::pair<double,double> > sssolmol = mutualScatteringListSolMol[index1][index2];
  for(int i=0;i<sssolmol.size();i++){
    std::pair<double,double> p = sssolmol[i];
    std::cout<<p.first<<" "<<p.second<<"\n";
  }
  std::vector<std::pair<double,double> > ssmol = mutualScatteringListSol[index1][index2];
  std::cout<<"mutual scatter section "<<index1<<" "<<index2<<"molecule:\n";
  for(int i=0;i<ssmol.size();i++){
    std::pair<double,double> p = ssmol[i];
    std::cout<<p.first<<" "<<p.second<<"\n";
  }
}

void hydrationShellMinimal::printTotalOneMutualScatter(int index1,int index2){
  std::vector<std::pair<double,double> > sssol = mutualScatteringListSol[index1][index2];
  std::cout<<"mutual scatter section "<<index1+1<<" "<<index2+1<<"Molecule:\n";
  double solsum = 0.0;
  for(int i=0;i<sssol.size();i++){
    std::pair<double,double> p = sssol[i];
    solsum = solsum +p.second;
  }
  std::cout<<solsum<<"\n";
}


double  hydrationShellMinimal::returnTotalOneMutualScatter(int index1,int index2){
   std::vector<std::pair<double,double> > ssmol = selfScatteringListSolMol[index1];
   double solsum;
   if(ssmol.size()>0){
     solsum = ssmol[0].second;
   }else{
     solsum =0;
   }
  /*for(int i=0;i<ssmol.size();i++){
    std::pair<double,double> p = ssmol[i];
    solsum = solsum +p.second;
    }*/
  return solsum;
}


void hydrationShellMinimal::printTotalScatter(){
  std::cout<<"total scatter solvent:\n";
  for(int i=0;i<totalScatterSol.size();i++){
    std::pair<double,double> p = totalScatterSol[i];
    std::cout<<p.first<<" "<<p.second<<"\n";
  }
  std::cout<<"total scatter solvent and molecule:\n";
  for(int i=0;i<totalScatterSolMol.size();i++){
    std::pair<double,double> p = totalScatterSolMol[i];
    std::cout<<p.first<<" "<<p.second<<"\n";
  }
  std::cout<<"total scatter molecule:\n";
  for(int i=0;i<totalScatterMol.size();i++){
    std::pair<double,double> p = totalScatterMol[i];
    std::cout<<p.first<<" "<<p.second<<"\n";
  }
}

void hydrationShellMinimal::writeTotalScatterToFile(const char* filename){
  std::ofstream myfile;
  myfile.open(filename);
  if(myfile.is_open()){
    myfile<<"total scatter solvent:\n";
    for(int i=0;i<totalScatterSol.size();i++){
      std::pair<double,double> p = totalScatterSol[i];
       myfile<<p.first<<" "<<p.second<<"\n";
    }
    myfile<<"total scatter solvent and molecule:\n";
    for(int i=0;i<totalScatterSolMol.size();i++){
      std::pair<double,double> p = totalScatterSolMol[i];
      myfile<<p.first<<" "<<p.second<<"\n";
    }
    myfile<<"total scatter molecule:\n";
    for(int i=0;i<totalScatterMol.size();i++){
      std::pair<double,double> p = totalScatterMol[i];
      myfile<<p.first<<" "<<p.second<<"\n";
    }
  }else{
    std::cout<<"scatter file failed to open \n";
  }
}


std::vector<std::vector<std::vector<int> > > hydrationShellMinimal::getIndicies(int index){
  std::vector<std::vector<std::vector<int> > > il = indexList[index];
  return il;
}

std::vector<std::vector<int> > hydrationShellMinimal::getTTS(int index){
  std::vector<std::vector<int> > il = allTruthTables[index];
  return il;
}



void hydrationShellMinimal::printOnePod(int index){
  std::cout<<"for pod "<<index<<"\n";
  hydTubeCentre(index).printPoint();
  getDirec(index).printPoint();
  getHydTan(index).printPoint();
  getHydNorm(index).printPoint();
  getHydBinorm(index).printPoint();
  std::cout<<0.5*hydTubeLength(index)<<"\n";
  std::cout<<"solvents"<<"\n";
  std::cout<<"\n";
  std::vector<std::vector<point> > solvents = allSegments[index];
  for(int i=0;i<solvents.size();i++){
    std::cout<<"segment "<<i<<"\n";
    for(int j=0;j<solvents[0].size();j++){
      solvents[i][j].printPoint();
    }
  }
  std::cout<<"\n";
  std::cout<<"truth table\n";
  std::cout<<"\n";
  std::vector<std::vector<int> > tt = allTruthTables[index];
  for(int i=0;i<tt.size();i++){
    std::cout<<"segment "<<i<<"\n";
    for(int j=0;j<tt[0].size();j++){
      std::cout<<tt[i][j]<<"\n";
    }
  }
  std::cout<<"\n";
  std::cout<<"Overlap Indicies\n";
  std::cout<<"\n";
  std::cout<<"\n";
  std::vector<std::vector<std::vector<int> > > il = indexList[index];
  for(int i=0;i<il.size();i++){
    std::cout<<"segment "<<i<<"\n";
    for(int j=0;j<il[0].size();j++){
      for(int k=0;k<il[i][j].size();k++){
	std::cout<<il[i][j][k]<<" ";
      }
      std::cout<<"\n";
    }
  }
  std::cout<<"\n";
  std::cout<<"\n";
  std::cout<<"minimum distances \n";
  std::cout<<"\n";
  std::cout<<"\n";
  std::cout<<minimumDistances[index].second<<" \n";
  std::cout<<" list\n\n";
  /*for(int i=0;i<minimumDistances[index].first.size();i++){
   std::cout<<minimumDistances[index].first[i].first<<" "<<minimumDistances[index].first[i].second<<"\n";
   }*/
  // std::cout<<"\n";
  //std::cout<<"\n";
  //std::cout<<" max change is "<<maxDistChange<<"\n";
  //std::cout<<"\n";
  //std::cout<<"\n";
}

void hydrationShellMinimal::printShellFrame(int index){
  std::cout<<"for pod "<<index<<"\n";
  hydTubeCentre(index).printPoint();
  getDirec(index).printPoint();
  getHydTan(index).printPoint();
  getHydNorm(index).printPoint();
  getHydBinorm(index).printPoint();
  std::cout<<0.5*hydTubeLength(index)<<"\n";
  std::cout<<"solvents"<<"\n";
  std::cout<<"\n";
  std::vector<std::vector<point> > solvents = allSegments[index];
  for(int i=0;i<solvents.size();i++){
    std::cout<<"segment "<<i<<"\n";
    for(int j=0;j<solvents[0].size();j++){
      solvents[i][j].printPoint();
    }
  }
}

void hydrationShellMinimal::printIndicies(int index){
  std::cout<<"for pod "<<index<<"\n";
  std::cout<<"Overlap Indicies\n";
  std::cout<<"\n";
  std::vector<std::vector<std::vector<int> > > il = indexList[index];
  for(int i=0;i<il.size();i++){
    std::cout<<"segment "<<i<<"\n";
    for(int j=0;j<il[0].size();j++){
      for(int k=0;k<il[i][j].size();k++){
	std::cout<<il[i][j][k]<<" ";
      }
      std::cout<<"\n";
    }
  }
}

void hydrationShellMinimal::printMinDists(int index){
  std::cout<<"for pod "<<index<<"\n";
  std::cout<<"\n";
  std::cout<<"\n";
  std::cout<<"minimum distances \n";
  std::cout<<"\n";
  std::cout<<"\n";
  std::cout<<minimumDistances[index].second<<" \n";
  std::cout<<" list\n\n";
  for(int i=0;i<minimumDistances[index].first.size();i++){
    std::cout<<minimumDistances[index].first[i].first<<" "<<minimumDistances[index].first[i].second<<"\n";
  }
  std::cout<<"\n";
  std::cout<<"\n";
  std::cout<<" max change is "<<maxDistChange<<"\n";
  std::cout<<"\n";
  std::cout<<"\n";
}

void hydrationShellMinimal::printIndiciesPrev(int index){
  std::cout<<"for pod "<<index<<"\n";
  std::cout<<"Overlap Indicies\n";
  std::cout<<"\n";
  std::vector<std::vector<std::vector<int> > > il = indexListPrev[index];
  for(int i=0;i<il.size();i++){
    std::cout<<"segment "<<i<<"\n";
    for(int j=0;j<il[0].size();j++){
      for(int k=0;k<il[i][j].size();k++){
	std::cout<<il[i][j][k]<<" ";
      }
      std::cout<<"\n";
    }
  }
   std::cout<<"\n";
}

void hydrationShellMinimal::writeHydrationShellToFile(const char* filename){
  std::ofstream ofl;
  ofl.open(filename);
  if(ofl.is_open()){
    for(int k=0;k<allSegments.size();k++){
      for(int i=0;i<allSegments[k].size();i++){
	for(int j=0;j<allSegments[k][0].size();j++){
	  if(allTruthTables[k][i][j]==1){
	    ofl<<allSegments[k][i][j].getX()<<" "<<allSegments[k][i][j].getY()<<" "<<allSegments[k][i][j].getZ()<<"\n";
	  }
	}
      }
    }
  }
}

int hydrationShellMinimal::getNoBins(double kmin,double kmax){
  std::vector<point> hypts;
  for(int k=0;k<allSegments.size();k++){
    for(int i=0;i<allSegments[k].size();i++){
      for(int j=0;j<allSegments[k][0].size();j++){
	if(allTruthTables[k][i][j]==1){
	  hypts.push_back(allSegments[k][i][j]);
	}
      }
    }
  }
  double max = 0.0;
  double dst;
  for(int i=0;i<hypts.size()-1;i++){
    for(int j=0;j<hypts.size();j++){
      dst = hypts[i].eDist(hypts[j]);
      if(dst>max){
	max =dst;
      }
    }
  }
  std::cout<<"twice max dist "<<2.0*max<<"\n";
  int noBins =  std::round(2.0*max*(kmax-kmin)*0.31830988618);
  return noBins;
}

int hydrationShellMinimal::getNoScatPoints(){
  return totalScatterMol.size();
}

void hydrationShellMinimal::reverseOverlaps(int &podIndex,int &changeIndex){
  if(podIndex<changeIndex){
    for(int i=0;i<indexList[podIndex].size();i++){
      for(int j=0;j<indexList[podIndex][0].size();j++){
        std::vector<int> v  = indexList[podIndex][i][j];
        v.erase(std::remove_if(v.begin(), v.end(),[changeIndex](const int& x) {return x >= changeIndex;}), v.end());
        //if we have emptied add the 1000000 entry
        if(v.size()==0){
          v.push_back(1000000);
        }
        indexList[podIndex][i][j]=v;
      }
    }
  }else{
    for(int i=0;i<indexList[podIndex].size();i++){
      for(int j=0;j<indexList[podIndex][0].size();j++){
        std::vector<int> v  = indexList[podIndex][i][j];
         v.erase(std::remove_if(v.begin(), v.end(),[changeIndex](const int& x) {return x < changeIndex;}), v.end());
        if(v.size()==0){
          v.push_back(1000000);
        }
        indexList[podIndex][i][j]=v;
      }
    }
  }
}

void hydrationShellMinimal::updateHydration(int &changeIndex){
   // pre clean
   indexListPrev= indexList;
   for(int i=0;i<molsize;i++){
     reverseOverlaps(i,changeIndex);
   }
   for(int i=0;i<changeIndex;i++){
     //we need to ascertain where the ac
   for(int j=changeIndex;j<molsize;j++){
      if(j==i+1){
	      overlapPodstructJoined(i,j);
      }else{
	      overlapPodstruct(i,j,i,j);
      }
    }
  }
   for(int i =0;i<indexList.size();i++){
      for(int k=0;k<indexList[i].size();k++){
	  for(int l=0;l<indexList[i][k].size();l++){
	    sort(indexList[i][k][l].begin(),indexList[i][k][l].end());
	  }
      }
   } 
  //check for changes
   for(int i=0;i<molsize;i++){
     if(equal(indexList[i].begin(),indexList[i].end(),indexListPrev[i].begin())==true){
       hydrationChangeList[i]=false;
     }else{
       hydrationChangeList[i]=true;
     }
   } 
}

void hydrationShellMinimal::reverseOverlapsNew(int &podIndex,int &changeIndex,int &changeIndexUp,std::vector<int>  &sublistAbove,std::vector<int> &sublistBelow){
  // first remove from list podIndex all potential overlaps above
  if(podIndex<changeIndex){
    for(int i=0;i<indexList[podIndex].size();i++){
      for(int j=0;j<indexList[podIndex][0].size();j++){
	std::vector<int> v  = indexList[podIndex][i][j];
	std::vector<int> vcop  = indexList[podIndex][i][j];
	/*std::cout<<"before \n\n";
	for(int k=0;k<v.size();k++){
	  std::cout<<v[k]<<" ";
	}
	std::cout<<"\n";*/
	int offset=0;
	for(int k=0;k<vcop.size();k++){
	  if(std::find(sublistAbove.begin(),sublistAbove.end(),vcop[k])!= sublistAbove.end()){
	    v.erase(v.begin()+k-offset);
	    offset++;
	  }
	}
	if(v.size()==0){
	  v.push_back(1000000);
	  allTruthTables[podIndex][i][j]=1;
	}
	indexList[podIndex][i][j]=v;
      }
    }
  }else if(podIndex>changeIndexUp){
    for(int i=0;i<indexList[podIndex].size();i++){
      for(int j=0;j<indexList[podIndex][0].size();j++){
	std::vector<int> v  = indexList[podIndex][i][j];
	std::vector<int> vcop  = indexList[podIndex][i][j];
	int offset=0;
	for(int k=0;k<vcop.size();k++){
	  if(std::find(sublistBelow.begin(),sublistBelow.end(),vcop[k])!= sublistBelow.end()){
	    v.erase(v.begin()+k-offset);
	    offset++;
	  }
	}
	if(v.size()==0){
	  v.push_back(1000000);
	  allTruthTables[podIndex][i][j]=1;
	}
	indexList[podIndex][i][j]=v;
      }
    }
  }else{
    for(int i=0;i<indexList[podIndex].size();i++){
      for(int j=0;j<indexList[podIndex][0].size();j++){
	std::vector<int> v;
	v.push_back(1000000);
	allTruthTables[podIndex][i][j]=1;
	indexList[podIndex][i][j]=v;
      }
    }
  }
  
}


void hydrationShellMinimal::updateHydrationNew(int &changeIndexIn){
  //work out the index boundats
  int changeIndex =0;
  for(int l=0;l<changeIndexIn-1;l++){
    if(mol.getType(l)=="Helix" || mol.getType(l)=="Strand"){
      changeIndex = changeIndex+1;
    }else{
      changeIndex = changeIndex+mol.getUnitNo(l);
    }
  }
  int changeIndexUp;
  if(mol.getType(changeIndexIn) =="loop" || mol.getType(changeIndexIn) =="Strand") {
    changeIndexUp = changeIndex + mol.getUnitNo(changeIndexIn)+  mol.getUnitNo(changeIndexIn+1);
  }else{
    changeIndexUp = changeIndex + mol.getUnitNo(changeIndexIn);
  }
  indexListPrev= indexList;
  maxDistChange = mol.getMaxDistChange();
  std::vector<int> cutOffIndicies(mol.molSize()-1,0);
  std::vector<std::vector<int> > overlapUpperLists(mol.molSize());
  std::vector<std::vector<int> > overlapLowerLists(mol.molSize());
  std::cout<<"indicies ? "<<changeIndex<<" "<<changeIndexUp<<"\n";
  for(int i=0;i<changeIndexUp;i++){
    // for each section, below the change index find the new potential opverlap point
    // by taking intou account the maximum change in molecule distance. First the start guess
    cutOffIndicies[i] = minimumDistances[i].second;
    bool foundNewLimit= false;
    if(cutOffIndicies[i]< minimumDistances[i].first.size()){
      while(cutOffIndicies[i]< minimumDistances[i].first.size()&&foundNewLimit==false){
	//subtract the change in distance;
	minimumDistances[i].first[cutOffIndicies[i]].first = minimumDistances[i].first[cutOffIndicies[i]].first-maxDistChange;
	if(minimumDistances[i].first[cutOffIndicies[i]].first > 0){
	  foundNewLimit = true;
	}else{
	   cutOffIndicies[i]++;
	}	 
      }
    }
    if(cutOffIndicies[i] >= minimumDistances[i].first.size()){
      cutOffIndicies[i] = minimumDistances[i].first.size()-1;
    }
    //now we have found the new overlap
    for(int jv=0;jv<=cutOffIndicies[i];jv++){
      if(minimumDistances[i].first[jv].second>=changeIndex){
	overlapUpperLists[i].push_back(minimumDistances[i].first[jv].second);
	overlapLowerLists[minimumDistances[i].first[jv].second].push_back(i);
      }
    }
  }
    // now delete potential index changes
    for(int i=0;i<molsize;i++){
      reverseOverlapsNew(i,changeIndex,changeIndexUp,overlapUpperLists[i],overlapLowerLists[i]);
    }
    // First check if any points below the change and an above have changed overlap
    for(int i=0;i<changeIndex;i++){
      for(int jv=0;jv<=cutOffIndicies[i];jv++){
       if(minimumDistances[i].first[jv].second>= changeIndex){
	 //delete the appropriate records
	 if(minimumDistances[i].first[jv].second == i+1){
	   overlapPodstructJoined(i,i+1);
	 }else{
	   overlapPodstruct(i,minimumDistances[i].first[jv].second,i,jv+i+1);
	 }
       }
     }
   }
   // now we have to allow for the fact that the internal overlap of the particular section  might have changed
   for(int i=changeIndex;i<=changeIndexUp;i++){
     for(int jv=i+1;jv<molsize;jv++){
       if(jv==i+1){
	 overlapPodstructJoined(i,i+1);
       }else{
	 overlapPodstruct(i,jv,i,jv);
       }
     }
   }
   for(int i=0;i<changeIndex;i++){
     for(int jv=changeIndex;jv<=changeIndexUp;jv++){
       if(jv==i+1){
	 overlapPodstructJoined(i,i+1);
       }else{
	 overlapPodstruct(i,jv,i,jv);
       }
     }
   }
   for(int i =0;i<indexList.size();i++){
      for(int k=0;k<indexList[i].size();k++){
	for(int l=0;l<indexList[i][k].size();l++){
	  sort(indexList[i][k][l].begin(),indexList[i][k][l].end());
	  indexList[i][k][l].erase(unique(indexList[i][k][l].begin(),indexList[i][k][l].end()),indexList[i][k][l].end());
	}
      }
   }
     //check for changes
   for(int i=0;i<molsize;i++){
     if(equal(indexList[i].begin(),indexList[i].end(),indexListPrev[i].begin())==true){
       hydrationChangeList[i]=false;
     }else{
       hydrationChangeList[i]=true;
     }
   }
   /*Finally we re-order the minimsum distance vectors */
   for(int i=0;i<molsize-1;i++){
    std::sort(minimumDistances[i].first.begin(),minimumDistances[i].first.end(),[](const std::pair<double,int> &left, const std::pair<double,int> &right) {
    return left.first < right.first;
});   
    int k = 0;
    bool outOfRange =false;
    while(k<minimumDistances[i].first.size() && outOfRange==false){
      if(minimumDistances[i].first[k].first>0.0){
	outOfRange=true;
	minimumDistances[i].second=k;
      }else{
	k++;
      }
    }
  }
}

   

void hydrationShellMinimal::updateFlatSolList(int &changeIndex){
  int noflsol = 0;
  for(int k=0;k<allSegments.size();k++){
    std::vector<point> spts;
    for(int i=0;i<allSegments[k].size();i++){
      for(int j=0;j<allSegments[k][0].size();j++){
	if(allTruthTables[k][i][j]==1){
	  spts.push_back(allSegments[k][i][j]);
	}
      }
    }
    solPtsFlat[k] = spts;
    //std::cout<<"in update segment "<<k<<" "<<spts.size()<<"\n";
    noflsol = noflsol +spts.size();
  }
  //std::cout<<" number of flat points in update "<<noflsol<<"\n";
}


std::vector<std::pair<double,double> >  hydrationShellMinimal::updateBinDistancesMoleculeTest(int &nscat,double &kmin,double &kmax,int &changeIndexIn){
  isTooClose=false;
  std::vector<double> distancesNew;
  std::vector<double> distancesOld;
  std::vector<double> subDistances;
  isTooClose=false;
   int changeIndex =0;
  for(int l=0;l<changeIndexIn-1;l++){
    if(mol.getType(l)=="Helix" || mol.getType(l)=="Strand"){
      changeIndex = changeIndex+1;
    }else{
      changeIndex = changeIndex+mol.getUnitNo(l);
    }
  }
  int changeIndexUp;
  if(mol.getType(changeIndexIn) =="loop" || mol.getType(changeIndexIn) =="Strand") {
    changeIndexUp = changeIndex + mol.getUnitNo(changeIndexIn)+  mol.getUnitNo(changeIndexIn+1);
  }else{
    changeIndexUp = changeIndex + mol.getUnitNo(changeIndexIn);
  }
  for(int i=0;i<=changeIndexUp;i++){
    for(int j=changeIndex;j<helixpts.size();j++){
      for(int k=0;k<helixpts[i].size();k++){
	for(int l=0;l<helixpts[j].size();l++){
	  if(i!=j){
	    double dist = helixpts[i][k].eDist(helixpts[j][l]);
	    subDistances.push_back(helixpts[i][k].eDist(helixpts[j][l]));
	  }
	}
      }
      if(j>i){
	distancesNew.insert(distancesNew.end(),subDistances.begin(),subDistances.end());
	distancesOld.insert(distancesOld.end(),molmolDistSum[i][j-i-1].begin(),molmolDistSum[i][j-i-1].end());
	molmolDistSum[i][j-i-1]=subDistances;
      }else if(j<i){
	distancesNew.insert(distancesNew.end(),subDistances.begin(),subDistances.end());
	distancesOld.insert(distancesOld.end(),molmolDistSum[j][i-j-1].begin(),molmolDistSum[j][i-j-1].end());
	molmolDistSum[j][i-j-1]=subDistances;
      }
      subDistances.clear();
    }
  }
  if(distancesNew.size()!=0 && distancesOld.size()!=0){
    std::sort(distancesNew.begin(),distancesNew.end());
    std::sort(distancesOld.begin(),distancesOld.end());
    double binsize = 0.2;
    int nbins = int(ceil(maxdist/0.2));
    std::vector<double> binlist;
    for(int i =0;i<=nbins;i++){
      binlist.push_back(binsize*i);
    }
    int k=0;int kprev=0;
    int dlenNew = distancesNew.size();
    int dlenOld = distancesOld.size();
    std::vector<std::pair<int,int> > binSizesNew;
    //bin off the new distances
    std::pair<int,int> pr;
    for(int i=0;i<nbins;i++){
      while(distancesNew[k]<=binlist[i+1] && k<dlenNew){
	k++;
      }
      if((k-kprev)>=1){
	int noInBin = k-kprev;
	pr.first = i;
	pr.second = noInBin;
   distBinsMol[i] = distBinsMol[i]+noInBin;
	binSizesNew.push_back(pr);
	kprev=k;
      }
    }
    std::vector<std::pair<int,int> > binSizesOld;
  //bin off the old distances
    k=0;kprev=0;
    for(int i=0;i<nbins;i++){
      while(distancesOld[k]<=binlist[i+1] && k<dlenOld){
	k++;
      }
      if((k-kprev)>=1){
	int noInBin = k-kprev;
   pr.first=i;
	pr.second = noInBin;
  distBinsMol[i] = distBinsMol[i]-noInBin;
  if(binlist[i+1]<mutualDistCutOff && distBinsMol[i]>0){
    isTooClose=true;
  }
	binSizesOld.push_back(pr);
	kprev = k;
      }
    }
    // check for overlap 
    int l=0;
    while(binlist[l+1]<mutualDistCutOff){
    if(distBinsMol[l]>0){
      isTooClose=true;
    }
    l++;
    }
    
    double kstep = (kmax-kmin)/(nscat-1);
    std::vector<std::pair<double,double> > scatAdjust;
    std::pair<double,double> spr;
    for(int i =0;i<nscat;i++){
      double scatVal =0.0;
      for(int j=0;j<binSizesNew.size();j++){
	scatVal = scatVal + phaseList[i][binSizesNew[j].first]*binSizesNew[j].second;
      }
      for(int j=0;j<binSizesOld.size();j++){
	scatVal = scatVal - phaseList[i][binSizesOld[j].first]*binSizesOld[j].second;
      }
      spr.first = totalScatterMol[i].first;
    spr.second = scatVal;
    scatAdjust.push_back(spr);
    }
    return scatAdjust;
  }else{
    std::vector<std::pair<double,double> > scatAdjust;
    return scatAdjust;
  }
}

/*std::vector<std::pair<double,double> > hydrationShellMinimal::updateTest(int &nscat,double &kmin,double &kmax,int &changeIndexIn){
  std::vector<double> distancesNew;
  std::vector<double> distancesOld;
  std::vector<double> subDistances;
  //ne distances from below thwe
   int changeIndex =0;
  for(int l=0;l<changeIndexIn-1;l++){
    if(mol.getType(l)=="Helix" || mol.getType(l)=="Strand"){
      changeIndex = changeIndex+1;
    }else{
      changeIndex = changeIndex+mol.getUnitNo(l);
    }
  }
  int changeIndexUp;
  if(mol.getType(changeIndexIn) =="loop" || mol.getType(changeIndexIn) =="Strand") {
    changeIndexUp = changeIndex + mol.getUnitNo(changeIndexIn)+  mol.getUnitNo(changeIndexIn+1);
  }else{
    changeIndexUp = changeIndex + mol.getUnitNo(changeIndexIn);
  }
  for(int i=0;i<=changeIndexUp;i++){
    for(int j=changeIndex;j<solPtsFlat.size();j++){
      for(int k=0;k<solPtsFlat[i].size();k++){
	for(int l=0;l<solPtsFlat[j].size();l++){
	  if(i!=j){
	    subDistances.push_back(solPtsFlat[i][k].eDist(solPtsFlat[j][l]));
	  }
	}
      }
      if(j>i){
	distancesNew.insert(distancesNew.end(),subDistances.begin(),subDistances.end());
	distancesOld.insert(distancesOld.end(),solsolDistSum[i][j-i-1].begin(),solsolDistSum[i][j-i-1].end());
	solsolDistSum[i][j-i-1]=subDistances;
      }else if(j<i){
	distancesNew.insert(distancesNew.end(),subDistances.begin(),subDistances.end());
	distancesOld.insert(distancesOld.end(),solsolDistSum[j][i-j-1].begin(),solsolDistSum[j][i-j-1].end());
	solsolDistSum[j][i-j-1]=subDistances;
      }
      subDistances.clear();
    }
  }
  // contributions within the below change index section which occur if their shells have changed its size 
  for(int i=0;i<changeIndex-1;i++){
    for(int j=i+1;j<changeIndex;j++){
      if(hydrationChangeList[i]==1 || hydrationChangeList[j]==1){
	for(int k=0;k<solPtsFlat[i].size();k++){
	  for(int l=0;l<solPtsFlat[j].size();l++){
	    if(i!=j){
	      subDistances.push_back(solPtsFlat[i][k].eDist(solPtsFlat[j][l]));
	    }
	  }
	}
	if(j>i){
	  distancesNew.insert(distancesNew.end(),subDistances.begin(),subDistances.end());
	  distancesOld.insert(distancesOld.end(),solsolDistSum[i][j-i-1].begin(),solsolDistSum[i][j-i-1].end());
	  solsolDistSum[i][j-i-1]=subDistances;
	}else if(j<i){
	  distancesNew.insert(distancesNew.end(),subDistances.begin(),subDistances.end());
	  distancesOld.insert(distancesOld.end(),solsolDistSum[j][i-j-1].begin(),solsolDistSum[j][i-j-1].end());
	  solsolDistSum[j][i-j-1]=subDistances;
	}
	subDistances.clear();
      }
    }
  }
// now contibutions for the changes section withi itself, again only if the number of sehll points have changed
  for(int i=changeIndexUp+1;i<solPtsFlat.size()-1;i++){
    for(int j=i+1;j<solPtsFlat.size();j++){
      if(hydrationChangeList[i]==1 || hydrationChangeList[j]==1){
	for(int k=0;k<solPtsFlat[i].size();k++){
	  for(int l=0;l<solPtsFlat[j].size();l++){
	    if(i!=j){
	      subDistances.push_back(solPtsFlat[i][k].eDist(solPtsFlat[j][l]));
	    }
	  }
	}
	if(j>i){
	  distancesNew.insert(distancesNew.end(),subDistances.begin(),subDistances.end());
	  distancesOld.insert(distancesOld.end(),solsolDistSum[i][j-i-1].begin(),solsolDistSum[i][j-i-1].end());
	  solsolDistSum[i][j-i-1]=subDistances;
	}else if(j<i){
	  distancesNew.insert(distancesNew.end(),subDistances.begin(),subDistances.end());
	  distancesOld.insert(distancesOld.end(),solsolDistSum[j][i-j-1].begin(),solsolDistSum[j][i-j-1].end());
	  solsolDistSum[j][i-j-1]=subDistances;
	}
	subDistances.clear();
      }
    }
  }
  if(distancesNew.size()!=0 && distancesOld.size()!=0){
    std::sort(distancesNew.begin(),distancesNew.end());
    std::sort(distancesOld.begin(),distancesOld.end());
    double binsize = 0.2;
    int nbins = int(ceil(maxdist/0.2));
    std::vector<double> binlist;
    for(int i =0;i<=nbins;i++){
      binlist.push_back(binsize*i);
    }
    int k=0;int kprev=0;
    int dlenNew = distancesNew.size();
    int dlenOld = distancesOld.size();
    std::vector<std::pair<int,int> > binSizesNew;
    //bin off the new distances
    std::pair<int,int> pr;
    for(int i=0;i<nbins;i++){
      while(distancesNew[k]<=binlist[i+1] && k<dlenNew){
	k++;
      }
      if((k-kprev)>=1){
	int noInBin = k-kprev;
	pr.first = i;
	pr.second = noInBin;
	binSizesNew.push_back(pr);
	kprev=k;
      }
    }
    std::vector<std::pair<int,int> > binSizesOld;
  //bin off the old distances
    k=0;kprev=0;
    for(int i=0;i<nbins;i++){
      while(distancesOld[k]<=binlist[i+1] && k<dlenOld){
	k++;
      }
      if((k-kprev)>=1){
	int noInBin = k-kprev;
	pr.first = i;
	pr.second = noInBin;
	binSizesOld.push_back(pr);
	kprev = k;
      }
    }
    double kstep = (kmax-kmin)/(nscat-1);
    std::vector<std::pair<double,double> > scatAdjust;
    std::pair<double,double> spr;
    for(int i =0;i<nscat;i++){
      double scatVal =0.0;
      for(int j=0;j<binSizesNew.size();j++){
	scatVal = scatVal + phaseList[i][binSizesNew[j].first]*binSizesNew[j].second;
      }
      for(int j=0;j<binSizesOld.size();j++){
	scatVal = scatVal - phaseList[i][binSizesOld[j].first]*binSizesOld[j].second;
      }
      spr.first = totalScatterSol[i].first;
    spr.second = scatVal;
    scatAdjust.push_back(spr);
    }
    return scatAdjust;
  }else{
    std::vector<std::pair<double,double> > scatAdjust;
    return scatAdjust;
  }
}
*/
/*std::vector<std::pair<double,double> > hydrationShellMinimal::updateMoleculeTest(int &nscat,double &kmin,double &kmax,int &changeIndexIn){
  std::vector<double> distancesNew;
  std::vector<double> distancesOld;
  std::vector<double> subDistances;
  int changeIndex=0;
  for(int l=0;l<changeIndexIn-1;l++){
    if(mol.getType(l)=="Helix" || mol.getType(l)=="Strand"){
      changeIndex = changeIndex+1;
    }else{
      changeIndex = changeIndex+mol.getUnitNo(l);
    }
  }
  int changeIndexUp;
  if(mol.getType(changeIndexIn) =="loop" || mol.getType(changeIndexIn) =="Strand") {
    changeIndexUp = changeIndex + mol.getUnitNo(changeIndexIn)+  mol.getUnitNo(changeIndexIn+1);
  }else{
    changeIndexUp = changeIndex + mol.getUnitNo(changeIndexIn);
  }
  for(int i=0;i<changeIndex;i++){
    for(int j=changeIndex;j<solPtsFlat.size();j++){
      for(int k=0;k<helixpts[i].size();k++){
	for(int l=0;l<solPtsFlat[j].size();l++){
	    subDistances.push_back(helixpts[i][k].eDist(solPtsFlat[j][l]));
	}
      }
      distancesNew.insert(distancesNew.end(),subDistances.begin(),subDistances.end());
      distancesOld.insert(distancesOld.end(),molsolDistSum[i][j-i-1].begin(),molsolDistSum[i][j-i-1].end());
      molsolDistSum[i][j-i-1]=subDistances;
      subDistances.clear();
    }
  }
  // possible change of the  below change mols with the below change solvents
  for(int i=0;i<changeIndex-1;i++){
    for(int j=i+1;j<=changeIndex;j++){
      if(hydrationChangeList[j]==1||j>=changeIndex){
	for(int k=0;k<helixpts[i].size();k++){
	  for(int l=0;l<solPtsFlat[j].size();l++){
	    if(i!=j){
	      subDistances.push_back(helixpts[i][k].eDist(solPtsFlat[j][l]));
	    }
	  }
	}
	  distancesNew.insert(distancesNew.end(),subDistances.begin(),subDistances.end());
	  distancesOld.insert(distancesOld.end(),molsolDistSum[i][j-i-1].begin(),molsolDistSum[i][j-i-1].end());
	  molsolDistSum[i][j-i-1]=subDistances;
	subDistances.clear();
      }
    }
  }
  for(int i=changeIndex;i<=changeIndexUp;i++){
    for(int j=i+1;j<solPtsFlat.size();j++){
      for(int k=0;k<helixpts[i].size();k++){
	for(int l=0;l<solPtsFlat[j].size();l++){
	    subDistances.push_back(helixpts[i][k].eDist(solPtsFlat[j][l]));
	}
      }
      distancesNew.insert(distancesNew.end(),subDistances.begin(),subDistances.end());
      distancesOld.insert(distancesOld.end(),molsolDistSum[i][j-i-1].begin(),molsolDistSum[i][j-i-1].end());
      molsolDistSum[i][j-i-1]=subDistances;
      subDistances.clear();
    } 
  }
  // possible change of the  above change mols with the above change solvents
  for(int i=changeIndexUp+1;i<helixpts.size()-1;i++){
    for(int j=i+1;j<solPtsFlat.size();j++){
      if(hydrationChangeList[j]==1){
	for(int k=0;k<helixpts[i].size();k++){
	  for(int l=0;l<solPtsFlat[j].size();l++){
	      subDistances.push_back(helixpts[i][k].eDist(solPtsFlat[j][l]));
	  }
	}
	  distancesNew.insert(distancesNew.end(),subDistances.begin(),subDistances.end());
	  distancesOld.insert(distancesOld.end(),molsolDistSum[i][j-i-1].begin(),molsolDistSum[i][j-i-1].end());
	  molsolDistSum[i][j-i-1]=subDistances;
	subDistances.clear();
      }
    }
  }
  for(int i=0;i<changeIndex;i++){
    for(int j=changeIndex;j<helixpts.size();j++){
      for(int k=0;k<solPtsFlat[i].size();k++){
	for(int l=0;l<helixpts[j].size();l++){
	    subDistances.push_back(solPtsFlat[i][k].eDist(helixpts[j][l]));
	}
      }
	distancesNew.insert(distancesNew.end(),subDistances.begin(),subDistances.end());
	distancesOld.insert(distancesOld.end(),solmolDistSum[i][j-i-1].begin(),solmolDistSum[i][j-i-1].end());
	solmolDistSum[i][j-i-1]=subDistances;
 subDistances.clear();
    }
  }
  for(int i=0;i<changeIndex-1;i++){
    if(hydrationChangeList[i]==1){
      for(int j=i+1;j<changeIndex;j++){
	for(int k=0;k<solPtsFlat[i].size();k++){
	  for(int l=0;l<helixpts[j].size();l++){
	    subDistances.push_back(solPtsFlat[i][k].eDist(helixpts[j][l]));
	  }
	}
	distancesNew.insert(distancesNew.end(),subDistances.begin(),subDistances.end());
	distancesOld.insert(distancesOld.end(),solmolDistSum[i][j-i-1].begin(),solmolDistSum[i][j-i-1].end());
	solmolDistSum[i][j-i-1]=subDistances;
	subDistances.clear();
      }
    }
  }
  for(int i =changeIndex;i<=changeIndexUp;i++){
    for(int j=i+1;j<helixpts.size();j++){
      for(int k=0;k<solPtsFlat[i].size();k++){
	for(int l=0;l<helixpts[j].size();l++){
	    subDistances.push_back(solPtsFlat[i][k].eDist(helixpts[j][l]));
	}
      }
	distancesNew.insert(distancesNew.end(),subDistances.begin(),subDistances.end());
	distancesOld.insert(distancesOld.end(),solmolDistSum[i][j-i-1].begin(),solmolDistSum[i][j-i-1].end());
	solmolDistSum[i][j-i-1]=subDistances;
 subDistances.clear();
    }
  }
  for(int i=changeIndexUp+1;i<solPtsFlat.size()-1;i++){
    if(hydrationChangeList[i]==1){
      for(int j=i+1;j<helixpts.size();j++){
	for(int k=0;k<solPtsFlat[i].size();k++){
	  for(int l=0;l<helixpts[j].size();l++){
	      subDistances.push_back(solPtsFlat[i][k].eDist(helixpts[j][l]));
	  }
	}
	distancesNew.insert(distancesNew.end(),subDistances.begin(),subDistances.end());
	distancesOld.insert(distancesOld.end(),solmolDistSum[i][j-i-1].begin(),solmolDistSum[i][j-i-1].end());
	solmolDistSum[i][j-i-1]=subDistances;
	subDistances.clear();
      }
    }
  }
  if(distancesNew.size()!=0 && distancesOld.size()!=0){
    std::sort(distancesNew.begin(),distancesNew.end());
    std::sort(distancesOld.begin(),distancesOld.end());
    double binsize = 0.2;
    int nbins = int(ceil(maxdist/0.2));
    std::vector<double> binlist;
    for(int i =0;i<=nbins;i++){
      binlist.push_back(binsize*i);
    }
    int k=0;int kprev=0;
    int dlenNew = distancesNew.size();
    int dlenOld = distancesOld.size();
    std::vector<std::pair<int,int> > binSizesNew;
    //bin off the new distances
    std::pair<int,int> pr;
    for(int i=0;i<nbins;i++){
      while(distancesNew[k]<=binlist[i+1] && k<dlenNew){
	k++;
      }
      if((k-kprev)>=1){
	int noInBin = k-kprev;
      pr.first = i;
      pr.second = noInBin;
      binSizesNew.push_back(pr);
      kprev=k;
      }
    }
    std::vector<std::pair<int,int> > binSizesOld;
    //bin off the old distances
    k=0;kprev=0;
    for(int i=0;i<nbins;i++){
      while(distancesOld[k]<=binlist[i+1] && k<dlenOld){
	k++;
      }
      if((k-kprev)>=1){
	int noInBin = k-kprev;
	pr.first = i;
	pr.second = noInBin;
	binSizesOld.push_back(pr);
	kprev = k;
      }
    }
    double kstep = (kmax-kmin)/(nscat-1);
    std::vector<std::pair<double,double> > scatAdjust;
    std::pair<double,double> spr;
    
    for(int i =0;i<nscat;i++){
      double scatVal =0.0;
      for(int j=0;j<binSizesNew.size();j++){
	scatVal = scatVal + phaseList[i][binSizesNew[j].first]*binSizesNew[j].second;
      }
      for(int j=0;j<binSizesOld.size();j++){
	scatVal = scatVal - phaseList[i][binSizesOld[j].first]*binSizesOld[j].second;
      }
      spr.first = totalScatterSolMol[i].first;
      spr.second = scatVal;
      scatAdjust.push_back(spr);
    }
    return scatAdjust;
  }
  else{
    // here there was no need to update anything
    std::vector<std::pair<double,double> > scatAdjust;
    return scatAdjust;
  }
}
*/  

/*void hydrationShellMinimal::updateScatteringBinnedTest(double &kmin,double &kmax,int &changeIndex,int &changeIndexLow,int &changeIndexHigh){
  //set up
  int nscat = totalScatterSol.size();
  totalScatterSol = emptyScatterSol;
  totalScatterMol = emptyScatterMol;
  totalScatterSolMol = emptyScatterSolMol;
  if(std::abs(totalScatterSol[0].first - kmin)){
    for(int i=0;i<molsize;i++){
      if(hydrationChangeList[i]==true|| changeIndexLow<=i<=changeIndexHigh){
	selfScatteringSingleFastRT(nscat,kmin,kmax,i);
	addOneScatterSelf(nscat,i);
      }else{
	addOneScatterSelf(nscat,i);
      } 
    }
  }else{
    for(int i=0;i<molsize;i++){
      if(hydrationChangeList[i]==true|| changeIndexLow<=i<=changeIndexHigh){
	selfScatteringSingleFastRT(nscat,kmin,kmax,i);
	addOneScatterSelf(nscat,i);
      }else{
	addOneScatterSelf(nscat,i);
      } 
    }
  }
  //selfScatteringAllRT(nscat,kmin,kmax);
  //now the mutual scattering, first fill the distance bin
  std::vector<std::pair<double,double> > changeMutualScatMol = updateBinDistancesMoleculeTest(nscat,kmin,kmax,changeIndex);
  if(changeMutualScatMol.size()>0){
  addScatteringLists(mutualScatterMol,changeMutualScatMol);
  }
  std::vector<std::pair<double,double> > changeMutualScatSol = updateBinDistancesSolventTest(nscat,kmin,kmax,changeIndex);
  if(changeMutualScatSol.size()>0){
    addScatteringLists(mutualScatterSol,changeMutualScatSol);
  }
  std::vector<std::pair<double,double> > changeMutualScatSolMol = updateBinDistancesSolventMoleculeTest(nscat,kmin,kmax,changeIndex);
  if(changeMutualScatSolMol.size()>0){
    addScatteringLists(mutualScatterSolMol,changeMutualScatSolMol);
  }
  // now calculate the scattering
  addScatteringLists(totalScatterMol,mutualScatterMol);
  addScatteringLists(totalScatterSol,mutualScatterSol);
  addScatteringLists(totalScatterSolMol,mutualScatterSolMol);
  }*/

/*void hydrationShellMinimal::changeConfigurationBinned(int &changeIndex,double &kmin,double &kmax){
  int currHIndex =0;
  for(int l=0;l<changeIndex-1;l++){
    if(mol.getType(l)=="Helix" || mol.getType(l)=="Strand"){
      currHIndex = currHIndex+1;
    }else{
      currHIndex = currHIndex+mol.getUnitNo(l);
    }
  }
  int changeIndexUp;
  if(mol.getType(changeIndex) =="loop" || mol.getType(changeIndex) =="Strand") {
    changeIndexUp = currHIndex + mol.getUnitNo(changeIndex)+  mol.getUnitNo(changeIndex+1);
  }else{
    changeIndexUp = currHIndex + mol.getUnitNo(changeIndex);
  }
  mol.changeMoleculeSingle(changeIndex);
  updateTubeParamList(changeIndex,currHIndex);
  updateState(currHIndex);
  updateHydrationNew(changeIndex);
  updateHelices(changeIndex);
  updateFlatSolList(changeIndex);
  updateScatteringBinnedTest(kmin,kmax,changeIndex,currHIndex,changeIndexUp);
  }*/


void hydrationShellMinimal::changeConfigurationBinnedFull(int &changeIndex){
  mol.changeMoleculeSingle(changeIndex);
  getScatteringBinnedRTUpdate();
}

void hydrationShellMinimal::changeConfigurationBinnedFullFast(int &changeIndex){
  mol.changeMoleculeSingle(changeIndex);
  getScatteringBinnedRTUpdate(changeIndex);
}

void hydrationShellMinimal::changeConfigurationLocalBinnedFull(std::vector<int> &changeIndicies,double &variationSize){
  mol.changeMoleculeLocalSet(changeIndicies,variationSize);
  getScatteringBinnedRTUpdate();
}


void hydrationShellMinimal::getUpdate(int &changeIndex){
  mol.changeMoleculeSingle(changeIndex);
}

void hydrationShellMinimal::changeMoleculeSingle(int &index){
  mol.changeMoleculeSingle(index);
}

void hydrationShellMinimal::changeMoleculeSet(std::vector<int> &indicies){
  mol.changeMoleculeSet(indicies);
}


void hydrationShellMinimal::changeMoleculeSingleMulti(int &index,int secIn){
  mol.changeMoleculeSingleMulti(index,secIn);
}

void hydrationShellMinimal::changeMoleculeSetMulti(std::vector<int>  &indicies,int secIn){
  mol.changeMoleculeSetMulti(indicies,secIn);
}

void hydrationShellMinimal::changeMoleculeMultiRotate(double &angle,point &k,int secIn,point &transVec){
  mol.changeMoleculeMultiRotate(angle,k,secIn,transVec);
}

void hydrationShellMinimal::replicateMolecule(int &noReplications){
  mol.replicateMolecule(noReplications);
}

void hydrationShellMinimal::changeMoleculeLocal(int &index,double variationSize){
  mol.changeMoleculeLocal(index,variationSize);
}

void hydrationShellMinimal::changeMoleculeLocalSet(std::vector<int> &indicies,double variationSize){
  mol.changeMoleculeLocalSet(indicies,variationSize);
}

void hydrationShellMinimal::changeMoleculeSingle(int &index,std::vector<std::vector<point> > &cdsIn,std::vector<std::pair<std::string,int> > &nameSizeSubList){
  mol.changeMoleculeSingle(index,cdsIn,nameSizeSubList);
}

int hydrationShellMinimal::getSubsecSize(int sec){
  return mol.getSubsecSize(sec);
}

std::vector<std::vector<point> > hydrationShellMinimal::getCoordinates(){
  return mol.getCoordinates();
};

double hydrationShellMinimal::maxNeighbourDistSec(int &sec){
  return mol.maxNeighbourDistSec(sec);
}


std::vector<double>  hydrationShellMinimal::checkOverlapWithRad(double &wRad){
  return mol.checkOverlapWithRad(wRad);
}

void hydrationShellMinimal::getHydrophobicResidues(){
  mol.getHydrophobicResidues();
}

void hydrationShellMinimal::getCoiledCoilResidues(){
  mol.getCoiledCoilResidues();
}

double hydrationShellMinimal::getGlobalRadiusOfCurvature(){
  return mol.getGlobalRadiusOfCurvature();
}

/*double hydrationShellMinimal::getGlobalRadiusOfCurvatureBetweenSec(){
  return mol.getGlobalRadiusOfCurvatureBetweenSec();
}*/

double hydrationShellMinimal::coiledCoilPotential(){
  return mol.coiledCoilPotential();
}

double hydrationShellMinimal::coiledCoilPotentialBetween(){
  return mol.coiledCoilPotentialBetween();
}

double hydrationShellMinimal::getFitQualSbonds(std::vector<std::tuple<std::pair<int,int>,std::pair<int,int>,double> > &contactPairList,double weightCoeff){
  return mol.getFitQualSbonds(contactPairList,weightCoeff);
}


double hydrationShellMinimal::compareDistances(std::vector<std::vector<point> > &coords2){
  return mol.compareDistances(coords2);
}

std::string hydrationShellMinimal::getType(int &chainIndex,int &index){
  return mol.getType(chainIndex,index);
}

std::vector<std::pair<double,double> > hydrationShellMinimal::getSolScatter(){
  return totalScatterSol;
}

std::vector<std::pair<double,double> > hydrationShellMinimal::getMolScatter(){
  //return mutualScatterMol;
  return totalScatterMol;
}

std::vector<std::pair<double,double> > hydrationShellMinimal::getSolMolScatter(){
  return totalScatterSolMol;
}

std::vector<std::vector<std::pair<double,double> > > hydrationShellMinimal::getSelfSolScatter(){
  return selfScatteringListSol;
}

std::vector<std::vector<std::pair<double,double> > > hydrationShellMinimal::getSelfMolScatter(){
  return selfScatteringListMol;
}

std::vector<std::vector<std::pair<double,double> > > hydrationShellMinimal::getSelfSolMolScatter(){
  return selfScatteringListSolMol;
}

std::vector<std::vector<std::vector<std::pair<double,double> > > > hydrationShellMinimal::getMutualSolScatter(){
  return mutualScatteringListSol;
}

  std::vector<std::vector<std::vector<std::pair<double,double> > > > hydrationShellMinimal::getMutualMolScatter(){
  return mutualScatteringListMol;
}

  std::vector<std::vector<std::vector<std::pair<double,double> > > > hydrationShellMinimal::getMutualSolMolScatter(){
  return mutualScatteringListSolMol;
}


void hydrationShellMinimal::writeMoleculeToFile(const char* filename){
  mol.writeMoleculeToFile(filename);
}

void hydrationShellMinimal::writeHelixPointsToFile(const char* filename){
  std::ofstream ofile;
  ofile.open(filename);
  if(ofile.is_open()){
    for(int i=0;i<helixpts.size();i++){
      for(int j=0;j<helixpts[i].size();j++){
	ofile<<helixpts[i][j].getX()<<" "<<helixpts[i][j].getY()<<" "<<helixpts[i][j].getZ()<<"\n";
      }
    } 
  }else{
    std::cout<<"cannot write molecule to file";
  }
  ofile.close();
}

void hydrationShellMinimal::writeHydrationParametersaToFile(const char* filename){
  std::ofstream ofile;
  ofile.open(filename);
  if(ofile.is_open()){
    for(int i=0;i<direcList.size();i++){
      ofile<<midPointList[i].getX()<<" "<<midPointList[i].getY()<<" "<<midPointList[i].getZ()<<" "<<direcList[i].getX()<<" "<<direcList[i].getY()<<" "<<direcList[i].getZ()<<" "<<halfLengthList[i]<<"\n";
    }
  }else{
    std::cout<<"cannot write hydration parameters to file";
  }
  ofile.close();
}

void hydrationShellMinimal::getBackboneStats(){
  mol.getBackboneStats();
}

std::vector<std::pair<double,double> > hydrationShellMinimal::getKapTauVals(){
    return mol.getKapTauVals();
}
