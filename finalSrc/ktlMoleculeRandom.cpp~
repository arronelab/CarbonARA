

#include "ktlMoleculeRandom.h"

ktlMolecule::ktlMolecule(){
  kapvallink = 0.36932;
  kapvalbeta = 0.343782;
  kapvalalpha =0.380928;
  tauvallink = 0.13439;
  tauvalbeta = 0.123222;
  tauvalalpha = 0.145831;
  alvallink = 4.22128;
  alvalbeta = 4.13944;
  alvalalpha = 4.2329;
};

void ktlMolecule::setParams(double &rminIn,double &rmaxIn,double &lminIn){
 rmg.setParams(rminIn,rmaxIn,lminIn);
}

std::vector<int> ktlMolecule::getUnitNos(){
  return noPts;
}

std::vector<std::vector<point> > ktlMolecule::getTangents(){
  return tanlist;
};

std::vector<std::vector<point> > ktlMolecule::getNormals(){
  return normlist;
};

std::vector<std::vector<point> > ktlMolecule::getBinormals(){
  return binormlist;
};

std::vector<std::vector<point> > ktlMolecule::getCoordinates(){
  return coords;
};

std::vector<point> ktlMolecule::getCoordinatesSection(int i){
  return coords[i];
}

int ktlMolecule::getSubsecSize(int sec){
  // the minus 1 is becasue the labellingnumbers are 1,2,3 e.t.c but for chainList 1 is at index0
  return chainList[sec-1].second-chainList[sec-1].first+1;
}

 std::vector<std::vector<point> > ktlMolecule::getSubsecCoordinates(int &sec){
  std::vector<std::vector<point> >::const_iterator firstc =coords.begin()+chainList[sec].first;
  std::vector<std::vector<point> >::const_iterator secondc =coords.begin()+chainList[sec].second+1;
  std::vector<std::vector<point> > subcoords(firstc,secondc);
  return subcoords;
}

std::vector<std::pair<std::string,int> > ktlMolecule::getNameSizeListOfSection(int &sec){
  std::vector<std::pair<std::string,int> >::const_iterator first =nameSizeList.begin()+chainList[sec].first;
  std::vector<std::pair<std::string,int> >::const_iterator second =nameSizeList.begin()+chainList[sec].second+1;
  std::vector<std::pair<std::string,int> > subns(first,second);
  return subns;
}

double ktlMolecule::maxNeighbourDistSec(int &sec){
  // the minus 1 is becasue the labellingnumbers are 1,2,3 e.t.c but for chainList 1 is at index0
  double dmax=0.0;
  for(int i = chainList[sec-1].first;i<= chainList[sec-1].second;i++){
    for(int j=0;j<coords[i].size();j++){      double d;
      if(j==(coords[i].size()-1) && i<chainList[sec-1].second){
	d = coords[i][j].eDist(coords[i+1][0]);
      }else if(j==(coords[i].size()-1) && i==chainList[sec-1].second){
	d=0.0;
      }else{
	d = coords[i][j].eDist(coords[i][j+1]);
      }
      if(d>dmax){
	dmax=d;
	if(d>4.0){
	  //std::cout<<i<<" "<<j<<"\n";
	}
      }
    }
  }
  return dmax;
}

/*std::vector<double> ktlMolecule::getDistChanges(){
  return distChanges;
  }*/


double ktlMolecule::getCurvatureJoined(int index){
  double kapval;
  int sz =coords[index].size();
  if(sz<=2){
    kapval = kapvallink;
  }else if(3<sz<7){
    kapval = kapvalbeta;
  }else{
    kapval = kapvalalpha;
  }
  return kapval;
}


double ktlMolecule::getTorsionJoined(int index){
  double tauval;
  int sz =coords[index].size();
  if(sz<=2){
    tauval = tauvallink;
  }else if(3<sz<7){
    tauval = tauvalbeta;
  }else{
    tauval = tauvalalpha;
  }
  return tauval;
}


 int ktlMolecule::getUnitNo(int index){
   return coords[index].size();
}

point ktlMolecule::getTangent(int mindex,int subindex){
  return tanlist[mindex][subindex];
};

point ktlMolecule::getNormal(int mindex,int subindex){
  return normlist[mindex][subindex];
};

point ktlMolecule::getBinormal(int mindex,int subindex){
  return binormlist[mindex][subindex];
};

point ktlMolecule::getCoordinate(int mindex,int subindex){
  if(subindex==-1){
    return coords[mindex-1][coords[mindex-1].size()-1];
  }else{
    return coords[mindex][subindex];
  }
};


double ktlMolecule::getAlbeadJoined(int index){
  double length;
  int sz =coords[index].size();
  if(sz<=2){
    length = alvallink;
  }else if(3<sz<7){
    length = alvalbeta;
  }else{
    length = alvalalpha;
  }
  return length;
}

std::string ktlMolecule::getType(int &index){
  return nameSizeList[index].first;
}

std::string ktlMolecule::getType(int &chainNo,int &index){
  int fullIndex = chainList[chainNo-1].first;
  return nameSizeList[fullIndex+index].first;
}


/*double ktlMolecule::getDistChange(int index){
  return distChanges[index];
  }*/

double ktlMolecule::getMaxDistChange(){
  return maxDistChange;
}

/*double ktlMolecule::getMaxScatLength(){
  double maxLen =0.0;
  double len=0.0;
  for(int i=0;i<coords.size();i++){
    for(int j=i+1;j<coords.size();j++){
      len = coords[i].eDist(coords[j]);
      if(len>maxLen){
        maxLen = len;
      }
    }
  }
  return 2.0*maxLen;
  }*/

int ktlMolecule::noSecSize(){
  return coords.size();
}

int ktlMolecule::noChains(){
  return chainList.size();
}

int ktlMolecule::molSize(){
  int sz = 0;
  for(int i=0;i<coords.size();i++){
    if(nameSizeList[i].first=="Helix" || nameSizeList[i].first=="Strand"){
      sz = sz + 1;
    }else{
      sz = sz + coords[i].size();
    }
  }
  return sz;
}

int ktlMolecule::getNoAminos(){
  int sz = 0;
  for(int i=0;i<coords.size();i++){
    sz=sz+coords[i].size();
  }
  return sz;
}

std::pair<double,double> ktlMolecule::getMaxPossibleLength(){
  double maxL=0.0;
  double totLength=0.0;
  double dst;
  for(int i=0;i<coords.size()-1;i++){
    for(int j=i+1;j<coords.size();j++){
      dst =coords[i][0].eDist(coords[j][0]);
      if(j==i+1){
	totLength= totLength+dst;
      }
      if(dst>maxL){
	maxL = dst;
      }
    }
  }
  maxL = maxL*2.0;
  std::pair<double,double> mxpr;
  mxpr.first = totLength;
  mxpr.second = maxL;
  return mxpr;
}



void ktlMolecule::readInMol(const char* filename,int chainNo,int isRand){
  int npts;
  std::vector<point> testPolyHelix;
  char fieldloc[1000]={};
  strcpy(fieldloc,"/home/rdata/ktch24/c++Molecule/calphaData/");
  strcat(fieldloc,filename);
  strcat(fieldloc,"/");
  strcat(fieldloc,filename);
  strcat(fieldloc,"ch_");
  std::string Result;
  std::stringstream convert;
  convert <<chainNo;
  Result = convert.str();
  const char* nolabel = (char*)Result.c_str();
  strcat(fieldloc,nolabel);
  strcat(fieldloc,"nameList.dat");
  std::ifstream myfile;
  myfile.open(fieldloc);
  std::string output;
  double val,X,Y,Z;
  int power;
  int n =0;
  int nunits;
  int nosecs=0;
  if (myfile.is_open()) { 
    while(!myfile.eof()){
      n++;
      std::getline(myfile,output);
      std::istringstream iss(output);
      std::pair<std::string,double> stpr;
      int i =0;
      for(std::string output;iss>>output;){
	      if(i==0){
       stpr.first=output;
	    }else{
	      std::stringstream ss(output);
       	 ss>>nunits;
        	stpr.second = nunits;
	    }
	      i++;
      }
      if(stpr.first=="Strand" ||stpr.first=="Helix"||stpr.first=="loop"){
      nameSizeList.push_back(stpr);
      distChanges.push_back(0.0);
      }
    }
    myfile.close();  
  }else{
    std::cout<<"Curve data file failed to open\n";
  }
}

void ktlMolecule::readInMolGen(const char* filename){
  int npts;
  std::vector<point> testPolyHelix;
  std::ifstream myfile;
  myfile.open(filename);
  std::string output;
  double val,X,Y,Z;
  int power;
  int n =0;
  int nunits;
  int nosecs=0;
  if (myfile.is_open()) { 
    while(!myfile.eof()){
      n++;
      std::getline(myfile,output);
      std::istringstream iss(output);
      std::pair<std::string,double> stpr;
      int i =0;
      for(std::string output;iss>>output;){
	      if(i==0){
       stpr.first=output;
	    }else{
	      std::stringstream ss(output);
       	 ss>>nunits;
        	stpr.second = nunits;
	    }
	      i++;
      }
      if(stpr.first=="Strand" ||stpr.first=="Helix"||stpr.first=="loop"){
      nameSizeList.push_back(stpr);
      distChanges.push_back(0.0);
      }
    }
    myfile.close();  
    std::pair<int,int> p;
    p.first = 0;
    p.second = nameSizeList.size()-1;      
    chainList.push_back(p);
  }else{
    std::cout<<"Curve data file failed to open\n";
  }
}

void ktlMolecule::readInMolWithBackbone(const char* filename,int chainNo,int isRand){
  int npts;
  std::vector<point> testPolyHelix;
  char fieldloc[1000]={};
  strcpy(fieldloc,"/home/rdata/ktch24/c++Molecule/calphaData/");
  strcat(fieldloc,filename);
  strcat(fieldloc,"/");
  strcat(fieldloc,filename);
  strcat(fieldloc,"ch_");
  std::string Result;
  std::stringstream convert;
  convert <<chainNo;
  Result = convert.str();
  const char* nolabel = (char*)Result.c_str();
  strcat(fieldloc,nolabel);
  strcat(fieldloc,"nameList.dat");
  std::ifstream myfile;
  //std::cout<<fieldloc<<"\n";
  myfile.open(fieldloc);
  std::string output;
  double val,prevval,X,Y,Z;
  int power;
  int n =0;
  int nunits;
  int nosecs=0;
  std::vector<point> section;
  if (myfile.is_open()) { 
    while(!myfile.eof()){
      n++;
      std::getline(myfile,output);
      std::istringstream iss(output);
      iss>>val;
      if(val==0){
      std::pair<std::string,double> stpr;
      int i =0;
      std::istringstream iss2(output);
      for(std::string output;iss2>>output;){
	if(i==0){
	  stpr.first=output;
	}else{
	  std::stringstream ss(output);
	  ss>>nunits;
	  stpr.second = nunits;
	}
	i++;
      }
      if(stpr.first=="Strand" ||stpr.first=="Helix"||stpr.first=="loop"){
      nameSizeList.push_back(stpr);
       distChanges.push_back(0.0);
      }
      }else{
    	if(prevval!=val){
             point p(output);
             section.push_back(p);
         }else{
           if(section.size()>0){
             coords.push_back(section);
           }
           section.clear();
         }
      }
      prevval=val;  
    }
    
    myfile.close();  
  }else{
    std::cout<<"Curve data file failed to open\n";
  }
}

void ktlMolecule::readInMolWithSequenceLenJ(const char* filename,int chainNo,int isRand,const char* LenJFileName){
  int npts;
  char fieldloc[1000]={};
  strcpy(fieldloc,"/home/rdata/ktch24/c++Molecule/calphaData/");
  strcat(fieldloc,filename);
  strcat(fieldloc,"/");
  strcat(fieldloc,filename);
  strcat(fieldloc,"Sequence");
  std::string Result;
  std::stringstream convert;
  convert <<chainNo;
  Result = convert.str();
  const char* nolabel = (char*)Result.c_str();
  strcat(fieldloc,nolabel);
  strcat(fieldloc,".dat");
  std::ifstream myfile;
  myfile.open(fieldloc);
  std::string output;
  double val,prevval,X,Y,Z;
  int power;
  int n =0;
  bool isFirst=true;
  int nunits;
  int nosecs=0;
  std::string type;
  std::string prevType;
  if (myfile.is_open()) { 
    while(!myfile.eof()){
      std::getline(myfile,output);
     /* input should read as a letter (amino type, currently ignore) and then the type H (helix) S (strand) or C (loop)
      */
      if(output.size()>1){
	if(isFirst){
	  if(output.compare(8,1,"H")==0){
	    type = "Helix";
	  }else if(output.compare(8,1,"S")==0){
	    type = "Strand";
	  }else{
	    type = "loop";
	  }
	  prevType=type;
	  n++;
	  isFirst=false;
	}else{
	  if(output.compare(8,1,"H")==0){
	    type = "Helix";
	  }else if(output.compare(8,1,"S")==0){
	    type = "Strand";
	  }else{
	    type = "loop";
	  }
	  if(prevType==type){
	    // keep increasing the section size
	    n++;
	  }else{
	    //store the previous type and length
	    std::pair<std::string,int> stpr;
	    stpr.first = prevType;stpr.second = n;
	    nameSizeList.push_back(stpr);
	    distChanges.push_back(0.0);
	    prevType = type;
	    n=1;
	  }
	}
      }
    }
    myfile.close();
  }
     // add the last section
  std::pair<std::string,int> stpr;
  stpr.first = type;stpr.second = n;
  nameSizeList.push_back(stpr);
  distChanges.push_back(0.0);
  std::ifstream myfileLenJ;
  myfileLenJ.open(LenJFileName);
  //get past the first line
  std::string outputLenJ;
  std::getline(myfileLenJ,outputLenJ);
  std::vector<point> section; 
  if (myfileLenJ.is_open()){ 
    while(!myfileLenJ.eof()){
      std::getline(myfileLenJ,outputLenJ);
      if(outputLenJ.size()>1){
        if(outputLenJ.find("End chain ")==std::string::npos){
	  point p2(outputLenJ);
	  section.push_back(p2);
	}
      }else{
        if(section.size()>0){
	  coords.push_back(section);
	}
	section.clear();
      }
    }
    myfile.close();
  }else{
    std::cout<<"Curve data file failed to open\n";
  }
}


void ktlMolecule::readInMolWithBackboneLenJ(const char* filename,int chainNo,int isRand,const char* LenJFileName){
  int npts;
  std::vector<point> testPolyHelix;
  char fieldloc[1000]={};

  std::ifstream myfile;
  std::ifstream myfileLenJ;
  std::cout<<"man is bare in this routine "<<fieldloc<<"\n";
  myfile.open(fieldloc);
  myfileLenJ.open(LenJFileName);
  //get past the first line
  std::string output;
  std::string outputLenJ;
  std::getline(myfileLenJ,outputLenJ);
  double val,val2,prevval,X,Y,Z;
  int power;
  int n =0;
  int nunits;
  int nosecs=0;
  std::vector<point> section;
  if (myfile.is_open()){ 
    while(!myfile.eof()){
      n++;
      std::getline(myfile,output);
      std::istringstream iss(output);
      iss>>val;
      if(val==0){
	std::pair<std::string,int> stpr;
	int i =0;
	std::istringstream iss2(output);
	for(std::string output;iss2>>output;){
	  if(i==0){
	    stpr.first=output;
	  }else{
	    std::stringstream ss(output);
	    ss>>nunits;
	    stpr.second = nunits;
	  }
	  i++;
	}
	if(stpr.first=="Strand" ||stpr.first=="Helix"||stpr.first=="loop"){
	  nameSizeList.push_back(stpr);
	  distChanges.push_back(0.0);
	}
      }else{
    	if(prevval!=val){
	  //point p(output);
	     std::getline(myfileLenJ,outputLenJ);
	     point p2(outputLenJ);
             section.push_back(p2);
         }else{
           if(section.size()>0){
             coords.push_back(section);
	     std::getline(myfileLenJ,outputLenJ);
           }
           section.clear();
         }
      }
      prevval=val;  
    }
    myfile.close();
  }else{
    std::cout<<"Curve data file failed to open\n";
  }
}


/*void ktlMolecule::readInSequence(const char* filename,int chainNo){
  int npts;
  char fieldloc[1000]={};
  strcpy(fieldloc,"/home/rdata/ktch24/c++Molecule/calphaData/");
  strcat(fieldloc,filename);
  strcat(fieldloc,"/");
  strcat(fieldloc,filename);
  strcat(fieldloc,"Sequence");
  std::string Result;
  std::stringstream convert;
  convert <<chainNo;
  Result = convert.str();
  const char* nolabel = (char*)Result.c_str();
  strcat(fieldloc,nolabel);
  strcat(fieldloc,".dat");
  std::ifstream myfile;
  //std::cout<<"in sequence ? "<<fieldloc<<"\n";
  myfile.open(fieldloc);
  std::string output;
  double val,prevval,X,Y,Z;
  int power;
  int n =0;
  bool isFirst=true;
  int nunits;
  int nosecs=0;
  std::string type;
  std::string prevType;
  if (myfile.is_open()) { 
    while(!myfile.eof()){
      std::getline(myfile,output);
     // input should read as a letter (amino type, currently ignore) and then the type H (helix) S (strand) or C (loop)
      if(output.size()>1){
	if(isFirst){
	  if(output.compare(8,1,"H")==0){
	    type = "Helix";
	  }else if(output.compare(8,1,"S")==0){
	    type = "Strand";
	  }else{
	    type = "loop";
	  }
	  prevType=type;
	  n++;
	  isFirst=false;
	}else{
	  if(output.compare(8,1,"H")==0){
	    type = "Helix";
	  }else if(output.compare(8,1,"S")==0){
	    type = "Strand";
	  }else{
	    type = "loop";
	  }
	  if(prevType==type){
	    // keep increasing the section size
	    n++;
	  }else{
	    //store the previous type and length
	    std::pair<std::string,int> stpr;
	    stpr.first = prevType;stpr.second = n;
	    nameSizeList.push_back(stpr);
	    distChanges.push_back(0.0);
	    prevType = type;
	    n=1;
	  }
	}
      }
    }
    myfile.close();
    // add the last section
     std::pair<std::string,int> stpr;
     stpr.first = type;stpr.second = n;
     nameSizeList.push_back(stpr);
     distChanges.push_back(0.0);
  }else{
    std::cout<<"Curve data file failed to open\n";
  }
}
*/

void ktlMolecule::readInSequence(const char* filename,double &rmin,double &rmax,double &lmin){
  int npts;
  std::ifstream myfile;
  myfile.open(filename);
  std::string output;
  double val,prevval,X,Y,Z;
  int power;
  int n =0;
  bool isFirst=true;
  int nunits;
  int nosecs=0;
  std::string sequence;
  std::string predictions;
  std::string empty;
  std::string chainNo;
  std::string type,prevType;
  std::vector<std::string> aminoType;
  int noChains;
  if(myfile.is_open()) {
    // read in the first line which tells us howmany chains there are
    std::getline(myfile,chainNo);
    std::stringstream ss(chainNo);
    ss>>noChains;
    std::cout<<"number of chains "<<noChains<<"\n";
    distSetsSecs.resize(noChains);
    for(int i=1;i<=noChains;i++){
      std::pair<int,int> p;
      p.first=nameSizeList.size();
      // grab the amino sequence (currently not used)
      std::getline(myfile,empty);
      std::getline(myfile,sequence);
      //std::cout<<sequence<<"\n";
      // grab the predictions
      std::getline(myfile,empty);
      std::getline(myfile,predictions);
      //std::cout<<predictions<<"\n";
      // input should read as a letter (amino type, currently ignore) and then the type H (helix) S (strand) or C (loop)
      n=0;
      for(int i=0;i<predictions.size();i++){
	if(i==0){
	  if(predictions.compare(i,1,"-")==0){
	    type = "loop";
	  }else if(predictions.compare(i,1,"E")==0 || predictions.compare(i,1,"S")==0){
	    type = "Strand";
	  }else{
	    type = "Helix";
	  }
	  n=1;
	  prevType=type;
	  aminoType.push_back(sequence.substr(i,1));
	}else{
	  if(predictions.compare(i,1,"-")==0){
	    type = "loop";
	  }else if(predictions.compare(i,1,"E")==0 || predictions.compare(i,1,"S")==0){
	    type = "Strand";
	  }else{
	    type = "Helix";
	  }
	  if(prevType==type){
	    // keep increasing the section size
	    n++;
	    aminoType.push_back(sequence.substr(i,1));
	  }else{
	    //store the previous type and length, if it is a helix check for breaks
	    if(prevType=="Helix"){
	      // helix check for glycine or prolene
	      int currLower=0;
	      int currSize=0;
	      std::vector<std::string> aminoTypeSub;
	      for(int j=0;j<aminoType.size();j++){
		aminoTypeSub.push_back(aminoType[j]);
		currSize++;
		//helix breaker, currently removed
		//if ((((aminoType[j]=="G" || aminoType[j]=="P")&& (j>2 && j<(aminoType.size()-3))) &&currSize>3) || (j==aminoType.size()-1) ){
		if(j==aminoType.size()-1){
		  std::pair<std::string,int> stpr;
		  stpr.first = prevType;stpr.second = j-currLower+1;
		  currLower=j+1;
		  nameSizeList.push_back(stpr);	
		  aminoList.push_back(aminoTypeSub);
		  distChanges.push_back(0.0);
		  aminoTypeSub.clear();
		  currSize=0;
		  }
	      }
	      prevType=type;
	      aminoType.clear();
	      aminoType.push_back(sequence.substr(i,1));
	      n=1;
	    }else{
	      //strand or linker, no need to separate
	      std::pair<std::string,int> stpr;
	      stpr.first = prevType;stpr.second = n;
	      nameSizeList.push_back(stpr);
	      distChanges.push_back(0.0);
	      prevType = type;
	      aminoList.push_back(aminoType);
	      aminoType.clear();
	      // empty the aminoVector
	      aminoType.push_back(sequence.substr(i,1));
	      n=1;
	    }
	  }
	}
      }
      std::pair<std::string,int> stpr;
      stpr.first = prevType;stpr.second = n;
      nameSizeList.push_back(stpr);
      distChanges.push_back(0.0);
      aminoList.push_back(aminoType);
      p.second = nameSizeList.size()-1;
      chainList.push_back(p);
      /*std::cout<<chainList.size()<<"\n";
      for(int iv=chainList[chainList.size()-2].second+1;iv<nameSizeList.size();iv++){
	std::cout<<nameSizeList[iv].first<<" "<<nameSizeList[iv].second<<"\n";
	}*/
    }
  //don't forget to set the random moelcule parameets
  rmg.setParams(rmin,rmax,lmin);
}else{
    // here no sequence secondary predictions provided, will just assume one large linker
    std::cout<<"no sequence provided for chain number "<<chainNo<<"\n";
  }
}

void ktlMolecule::readInCoordinates(const char* filename){
  std::ifstream coordFile(filename);
  std::string line;
  //std::cout<<"here?\n";
  if(coordFile.is_open()){
    std::vector<point> secondarySec;
    for(int iv=0;iv<nameSizeList.size();iv++){
      //keep reading until we get a coordinate
      bool waitToFill=true;
      int noInSection = nameSizeList[iv].second;
      int currNoMols=0;
      //std::cout<<"section number "<<iv<<" "<<nameSizeList[iv].first<<" "<<noInSection<<" "<<nameSizeList.size()<<"\n";
      while((coordFile.eof()==false)&&waitToFill==true){
	std::getline(coordFile,line);
	if(line.find("chain")==std::string::npos){
	  if(line.size()>1){
	    point p(line);
	    currNoMols++;
	    secondarySec.push_back(p);
	    if(currNoMols==noInSection){
	      // std::cout<<" this big pushed back"<<secondarySec.size()<<"\n";
	      coords.push_back(secondarySec);
	      secondarySec.clear();
	      waitToFill=false;
	      //std::cout<<"here?\n";
	    }
	  } 
	}
      }
      //std::cout<<coords.size()<<" "<<coords[iv].size()<<" "<<nameSizeList.size()<<" "<<nameSizeList[iv].second<<"\n";
    }
  }else{
    std::cout<<"there is no valid coordinate file supplied\n";
  }
}

void ktlMolecule::readInSequenceWBackbone(const char* filename,int chainNo,const char* backbonename){
  int npts;
  std::cout<<"all up in this routine\n";
  char fieldloc[1000]={};
  strcpy(fieldloc,"/home/rdata/ktch24/c++Molecule/calphaData/");
  strcat(fieldloc,filename);
  strcat(fieldloc,"/");
  strcat(fieldloc,filename);
  strcat(fieldloc,"Sequence");
  std::string Result;
  std::stringstream convert;
  convert <<chainNo;
  Result = convert.str();
  const char* nolabel = (char*)Result.c_str();
  strcat(fieldloc,nolabel);
  strcat(fieldloc,".dat");
  std::ifstream myfile;
  std::string output;
  std::string outputBackbone;
  //std::cout<<"in sequence ? "<<fieldloc<<"\n";
  myfile.open(fieldloc);
  std::ifstream myfileBackbone;
  myfileBackbone.open(backbonename);
  std::vector<point> section;
  //waste the first line which is not a cooordinate
  std::getline(myfileBackbone,outputBackbone);
  double val,prevval,X,Y,Z;
  int power;
  int n =0;
  bool isFirst=true;
  int nunits;
  int nosecs=0;
  std::string type;
  std::string prevType;
  if (myfile.is_open()) { 
    while(!myfile.eof()){
      std::getline(myfile,output);
     /* input should read as a letter (amino type, currently ignore) and then the type H (helix) S (strand) or C (loop)
      */
      if(output.size()>1){
	if(isFirst){
	  if(output.compare(8,1,"H")==0){
	    type = "Helix";
	  }else if(output.compare(8,1,"S")==0){
	    type = "Strand";
	  }else{
	    type = "loop";
	  }
	  prevType=type;
	  n++;
	  isFirst=false;
	}else{
	  if(output.compare(8,1,"H")==0){
	    type = "Helix";
	  }else if(output.compare(8,1,"S")==0){
	    type = "Strand";
	  }else{
	    type = "loop";
	  }
	  if(prevType==type){
	    // keep increasing the section size
	    n++;
	  }else{
	    //store the previous type and length
	    std::pair<std::string,int> stpr;
	    stpr.first = prevType;stpr.second = n;
	    nameSizeList.push_back(stpr);
	    distChanges.push_back(0.0);
	    /*Now add the molecules*/
	    //std::cout<<"here mol size "<<n<<"\n";
	    for(int k=0;k<n;k++){
	      std::getline(myfileBackbone,outputBackbone);
	      point p2(outputBackbone);  
	      //p2.printPoint();
	      section.push_back(p2);	 
	    }
	    //there is a space in the file, burn an extra line
	    std::getline(myfileBackbone,outputBackbone);
	    coords.push_back(section);
	    section.clear();
	    prevType = type;
	    n=1;
	  }
	}
      }
    }
    myfile.close();
    // add the last section
     std::pair<std::string,int> stpr;
     stpr.first = type;stpr.second = n;
      for(int k=0;k<n;k++){
	std::getline(myfileBackbone,outputBackbone);
	point p2(outputBackbone);  
	//p2.printPoint();
	section.push_back(p2);	 
      }
      coords.push_back(section);
      section.clear();
     nameSizeList.push_back(stpr);
     distChanges.push_back(0.0);
  }else{
    std::cout<<"Curve data file failed to open\n";
  }
}





// reads in just the secondary structure, not the coordinates themselevs, used for a from PDB secondary stcuture prediction

void ktlMolecule::readInMolNmer(const char* filename,std::vector<int> &molIndicies,int isRand,double &rmin,double &rmax,double &lmin){
  rmg.setParams(rmin,rmax,lmin);
  int noChains;
  char fieldloc[1000]={};
  std::string filen;
  strcpy(fieldloc,"/home/rdata/ktch24/c++Molecule/calphaData/");
  strcat(fieldloc,filename);
  strcat(fieldloc,"/");
  if(molIndicies[0] == -1){
    // read in all sections of the composite, find the no of chains 
    DIR *dir;
    struct dirent *ent;
    molIndicies.clear();
    int l=1;
    if((dir = opendir(fieldloc))!=NULL){
	while((ent = readdir(dir))!= NULL){
	  filen = ent->d_name;
	  std::size_t found = filen.find("KTL");
	  if(found!=std::string::npos){
	    molIndicies.push_back(l);
	    l++;
	  }
	}
    }else{
      std::cout<<"failed to find directory, molecule won't read in\n";
    }
  }
  noChains =  molIndicies.size();
  std::pair<int,int> p;
  p.first = 0;
  for(int i=0;i<noChains;i++){
    readInMol(filename,molIndicies[i],isRand);
    p.second = nameSizeList.size()-1;      
    chainList.push_back(p);
    p.first = nameSizeList.size();
  }
  // fill the vector used to track distance changes 
}

void ktlMolecule::readInMolNmerAndCoordinates(const char* filename,std::vector<int> &molIndicies,int isRand,double &rmin,double &rmax,double &lmin){
  rmg.setParams(rmin,rmax,lmin);
  int noChains;
  char fieldloc[1000]={};
  std::string filen;
  strcpy(fieldloc,"/home/rdata/ktch24/c++Molecule/calphaData/");
  strcat(fieldloc,filename);
  strcat(fieldloc,"/");
  if(molIndicies[0] == -1){
    // read in all sections of the composite, find the no of chains 
    DIR *dir;
    struct dirent *ent;
    molIndicies.clear();
    int l=1;
    if((dir = opendir(fieldloc))!=NULL){
	while((ent = readdir(dir))!= NULL){
	  filen = ent->d_name;
	  std::size_t found = filen.find("KTL");
	  if(found!=std::string::npos){
	    molIndicies.push_back(l);
	    l++;
	  }
	}
    }else{
      std::cout<<"failed to find directory, molecule won't read in\n";
    }
  }
  noChains =  molIndicies.size();
  std::pair<int,int> p;
  p.first = 0;
  for(int i=0;i<noChains;i++){
    readInMolWithBackbone(filename,molIndicies[i],isRand);
    p.second = nameSizeList.size()-1;      
    chainList.push_back(p);
    p.first = nameSizeList.size();
  }
  // fill the vector used to track distance changes 
}

/*
void ktlMolecule::readInMolNmerSequence(const char* filename,std::vector<int> &molIndicies,double &rmin,double &rmax,double &lmin){
  rmg.setParams(rmin,rmax,lmin);
  int noChains;
  char fieldloc[1000]={};
  std::string filen;
  strcpy(fieldloc,"/home/rdata/ktch24/c++Molecule/calphaData/");
  strcat(fieldloc,filename);
  strcat(fieldloc,"/");
  if(molIndicies[0] == -1){
    // read in all sections of the composite, find the no of chains 
    DIR *dir;
    struct dirent *ent;
    molIndicies.clear();
    int l=1;
    if((dir = opendir(fieldloc))!=NULL){
	while((ent = readdir(dir))!= NULL){
	  filen = ent->d_name;
	  std::size_t found = filen.find("KTL");
	  if(found!=std::string::npos){
	    molIndicies.push_back(l);
	    l++;
	  }
	}
    }else{
      std::cout<<"failed to find directory, molecule won't read in\n";
    }
  }
  noChains =  molIndicies.size();
  std::pair<int,int> p;
  p.first = 0;
  for(int i=0;i<noChains;i++){
    readInSequence(filename,molIndicies[i]);
    p.second = nameSizeList.size()-1;      
    chainList.push_back(p);
    p.first = nameSizeList.size();
  }
}
*/

void ktlMolecule::readInMolNmerSequenceWBackbone(const char* filename,const char* backboneName,std::vector<int> &molIndicies,double &rmin,double &rmax,double &lmin){
  rmg.setParams(rmin,rmax,lmin);
  int noChains;
  char fieldloc[1000]={};
  std::string filen;
  strcpy(fieldloc,"/home/rdata/ktch24/c++Molecule/calphaData/");
  strcat(fieldloc,filename);
  strcat(fieldloc,"/");
  if(molIndicies[0] == -1){
    // read in all sections of the composite, find the no of chains 
    DIR *dir;
    struct dirent *ent;
    molIndicies.clear();
    int l=1;
    if((dir = opendir(fieldloc))!=NULL){
	while((ent = readdir(dir))!= NULL){
	  filen = ent->d_name;
	  std::size_t found = filen.find("KTL");
	  if(found!=std::string::npos){
	    molIndicies.push_back(l);
	    l++;
	  }
	}
    }else{
      std::cout<<"failed to find directory, molecule won't read in\n";
    }
  }
  noChains =  molIndicies.size();
  std::pair<int,int> p;
  p.first = 0;
  for(int i=0;i<noChains;i++){
    readInSequenceWBackbone(filename,molIndicies[i],backboneName);
    p.second = nameSizeList.size()-1;      
    chainList.push_back(p);
    p.first = nameSizeList.size();
  }
  // fill the vector used to track distance changes 
}



void ktlMolecule::readInMolNmerWithBackbone(const char* filename,std::vector<int> &molIndicies,int isRand,double &rmin,double &rmax,double &lmin){
  rmg.setParams(rmin,rmax,lmin);
  int noChains;
  char fieldloc[1000]={};
  std::string filen;
  strcpy(fieldloc,"/home/rdata/ktch24/c++Molecule/calphaData/");
  strcat(fieldloc,filename);
  strcat(fieldloc,"/");
  if(molIndicies[0] == -1){
    // read in all sections of the composite, find the no of chains 
    DIR *dir;
    struct dirent *ent;
    molIndicies.clear();
    int l=1;
    if((dir = opendir(fieldloc))!=NULL){
      while((ent = readdir(dir))!= NULL){
	  filen = ent->d_name;
	  std::size_t found = filen.find("nameList");
	  if(found!=std::string::npos){
	    molIndicies.push_back(l);
	    l++;
	  }
	}
    }else{
      std::cout<<"failed to find directory, molecule won't read in\n";
    }
  }
  noChains =  molIndicies.size();
  std::pair<int,int> p;
  p.first = 0;
  for(int i=0;i<noChains;i++){
    readInMolWithBackbone(filename,molIndicies[i],isRand);
    p.second = nameSizeList.size()-1;      
    chainList.push_back(p);
    p.first = nameSizeList.size();
  }
  // fill the vector used to track distance changes
   std::vector<int> overLapList = checkOverlap(coords);
  std::vector<int> currOverLapList= overLapList;
  std::random_device rdev{};
  std::default_random_engine generator1{rdev()};
  int l=0;
  std::vector<std::pair<std::string,int> > nameSizeListTemp= nameSizeList;
  while(l<20000 && currOverLapList.size()>0){
    std::vector<std::vector<point> > tempCoords = coords;
    if(currOverLapList.size()>1){
      std::uniform_int_distribution<int> distributionR(0,currOverLapList.size()-1);
      int iv= distributionR(generator1);
      if(currOverLapList[iv]==0){
	changeMoleculeSingle(currOverLapList[1],tempCoords,nameSizeListTemp);  
      }else{
	changeMoleculeSingle(currOverLapList[iv],tempCoords,nameSizeListTemp);  
      }
    }else{
      changeMoleculeSingle(currOverLapList[0],tempCoords,nameSizeListTemp);
    }
    overLapList = checkOverlap(tempCoords);
    if(overLapList.size()<currOverLapList.size()){
      currOverLapList=overLapList;
      coords = tempCoords;
    }
    l++;
  }
}


void ktlMolecule::readInMolNmerWithBackboneLenJ(const char* filename,std::vector<int> &molIndicies,int isRand,double &rmin,double &rmax,double &lmin,const char* filenameLenJ){
  rmg.setParams(rmin,rmax,lmin);
  int noChains;
  char fieldloc[1000]={};
  std::string filen;
  strcpy(fieldloc,"/home/rdata/ktch24/c++Molecule/calphaData/");
  strcat(fieldloc,filename);
  strcat(fieldloc,"/");
  if(molIndicies[0] == -1){
    // read in all sections of the composite, find the no of chains 
    DIR *dir;
    struct dirent *ent;
    molIndicies.clear();
    int l=1;
    if((dir = opendir(fieldloc))!=NULL){
      while((ent = readdir(dir))!= NULL){
	  filen = ent->d_name;
	  std::size_t found = filen.find("nameList");
	  if(found!=std::string::npos){
	    molIndicies.push_back(l);
	    l++;
	  }
	}
    }else{
      std::cout<<"failed to find directory, molecule won't read in\n";
    }
  }
  noChains =  molIndicies.size();
  std::pair<int,int> p;
  p.first = 0;
  for(int i=0;i<noChains;i++){
    readInMolWithBackboneLenJ(filename,molIndicies[i],isRand,filenameLenJ);
    p.second = nameSizeList.size()-1;      
    chainList.push_back(p);
    p.first = nameSizeList.size();
  }
}

void ktlMolecule::readInMolNmerWithSequenceLenJ(const char* filename,std::vector<int> &molIndicies,int isRand,double &rmin,double &rmax,double &lmin,const char* filenameLenJ){
  rmg.setParams(rmin,rmax,lmin);
  int noChains;
  char fieldloc[1000]={};
  std::string filen;
  strcpy(fieldloc,"/home/rdata/ktch24/c++Molecule/calphaData/");
  strcat(fieldloc,filename);
  strcat(fieldloc,"/");
  if(molIndicies[0] == -1){
    // read in all sections of the composite, find the no of chains 
    DIR *dir;
    struct dirent *ent;
    molIndicies.clear();
    int l=1;
    if((dir = opendir(fieldloc))!=NULL){
      while((ent = readdir(dir))!= NULL){
	  filen = ent->d_name;
	  std::size_t found = filen.find("nameList");
	  if(found!=std::string::npos){
	    molIndicies.push_back(l);
	    l++;
	  }
	}
    }else{
      std::cout<<"failed to find directory, molecule won't read in\n";
    }
  }
  noChains =  molIndicies.size();
  std::pair<int,int> p;
  p.first = 0;
  for(int i=0;i<noChains;i++){
    readInMolWithSequenceLenJ(filename,molIndicies[i],isRand,filenameLenJ);
    p.second = nameSizeList.size()-1;      
    chainList.push_back(p);
    p.first = nameSizeList.size();
  }
}

void ktlMolecule::getHydrophobicResidues(){
  for(int i=0;i<aminoList.size();i++){
    for(int j=0;j<aminoList[i].size();j++){
      if(aminoList[i][j]=="A"||aminoList[i][j]=="I"||aminoList[i][j]=="L"||aminoList[i][j]=="M"||aminoList[i][j]=="F"||aminoList[i][j]=="V"||aminoList[i][j]=="P"||aminoList[i][j]=="G"){
	std::pair<int,int> pr;
	pr.first = i;pr.second = j;
	hydroPhobicList.push_back(pr);
      }
    }
  }
}


std::vector<double> ktlMolecule::getHydrophobicDistance(std::vector<std::vector<point> > &solventList,double &maxSolDist){
  // maxSolDist is the size of the spehere I check for solvents in 
  std::vector<double> meanHydroAminoDists;
  for(int i=0;i<hydroPhobicList.size();i++){
    std::pair<int,int> pr = hydroPhobicList[i];
    point coord = coords[pr.first][pr.second];
    double solDist = 0.0;
    int noDistances=0;
    // loop over all solvents, I don't care where they come from.
    for(int k=0;k<solventList.size();k++){
      for(int l=0;l<solventList[k].size();l++){
	double sddist = coord.eDist(solventList[k][l]);
	if(sddist<maxSolDist){
	  solDist = solDist + sddist;
	  noDistances++;
	}
      }
    }
    if(noDistances>0){
      solDist = solDist/double(noDistances);
    }
    //std::cout<<solDist<<" "<<noDistances<<"\n";
    meanHydroAminoDists.push_back(double(noDistances)/double(hydroPhobicList.size()));
  }
  return meanHydroAminoDists;
}


void ktlMolecule::getCoiledCoilResidues(){
  for(int i=0;i<aminoList.size();i++){
    for(int j=0;j<aminoList[i].size();j++){
      if(aminoList[i][j]=="L"||aminoList[i][j]=="I"||aminoList[i][j]=="G"){
	//check if the amino belongs to a helix
	if(nameSizeList[i].first=="Helix"){
	  std::pair<int,int> pr;
	  pr.first = i;pr.second = j;
	  coiledCoilList.push_back(pr);
	}
      }
    }
  }
  /* std::cout<<"in coiled coild routine "<<coiledCoilList.size()<<"\n";
  for(int i=0;i<coiledCoilList.size();i++){
    std::pair<int,int> pr = coiledCoilList[i];
    std::cout<<"in section "<<pr.first<<" of type "<<nameSizeList[pr.first].first<<" at residue "<<pr.second<<"\n";
    }*/
}

double ktlMolecule::coiledCoilPotential(){
  double potentialOut=0.0;
  int noCalcs=0;
  for(int i=0;i<coiledCoilList.size()-1;i++){
    for(int j=i+1;j<coiledCoilList.size();j++){
      if(coiledCoilList[i].first!=coiledCoilList[j].first){
	point p1=coords[coiledCoilList[i].first][coiledCoilList[i].second];
	point p2=coords[coiledCoilList[j].first][coiledCoilList[j].second];
	double d = p1.eDist(p2);d=d-7.3;
  double dsix = d*d;
	potentialOut = potentialOut + dsix;
       noCalcs++;
      }
    }
  }
  if(noCalcs>1){
  potentialOut= potentialOut/double(noCalcs);
  }
  return potentialOut;
}


double ktlMolecule::coiledCoilPotentialBetween(int &secNo){
  double potentialOut=0.0;
  int noCalcs=0;
  int firstIndex =chainList[secNo].first;
  int secondIndex = chainList[secNo].second;
  // get the sublist of hydrophobic residues in this section 
  std::vector<std::pair<int,int> > coiledCoilSubList;
  for(int i=0;i<coiledCoilList.size();i++){
      if((coiledCoilList[i].first>=firstIndex)&&(coiledCoilList[i].first<=secondIndex)){
        coiledCoilSubList.push_back(coiledCoilList[i]);
      }
  }
  std::vector<std::pair<int,int> > coiledCoilNotSubList;
  for(int i=0;i<coiledCoilList.size();i++){
      if((coiledCoilList[i].first<firstIndex)||(coiledCoilList[i].first>secondIndex)){
        coiledCoilNotSubList.push_back(coiledCoilList[i]);
      }
  }
  for(int i=0;i<coiledCoilSubList.size();i++){
    for(int j=0;j<coiledCoilNotSubList.size();j++){
      point p1=coords[coiledCoilList[i].first][coiledCoilSubList[i].second];
      point p2=coords[coiledCoilNotSubList[j].first][coiledCoilNotSubList[j].second];
      double d = p1.eDist(p2);d=d-7.3;
      double dsix = d*d;
      potentialOut = potentialOut + dsix;
      noCalcs++;
    }
  }
  if(noCalcs>1){
    potentialOut= potentialOut/double(noCalcs);
  }
  return potentialOut;
}

double ktlMolecule::coiledCoilPotentialBetween(){
  double coiledCoilPotential=0.0;
  int noChains = chainList.size();
  //std::cout<<"no chains are ? "<<noChains<<"\n";
  for(int i=0;i<noChains;i++){
    coiledCoilPotential = coiledCoilPotential + 0.5*coiledCoilPotentialBetween(i);
  }
  return coiledCoilPotential;
}


void ktlMolecule::getPolarResidues(){
  for(int i=0;i<aminoList.size();i++){
    for(int j=0;j<aminoList[i].size();j++){
      if(aminoList[i][j]=="Q"||aminoList[i][j]=="N"||aminoList[i][j]=="H"||aminoList[i][j]=="S"||aminoList[i][j]=="T"||aminoList[i][j]=="Y"||aminoList[i][j]=="C"){
	std::pair<int,int> pr;
	pr.first = i;pr.second = j;
	polarList.push_back(pr);
      }
    }
  }
}

void ktlMolecule::getPositiveResidues(){
  for(int i=0;i<aminoList.size();i++){
    for(int j=0;j<aminoList[i].size();j++){
      if(aminoList[i][j]=="K"||aminoList[i][j]=="Y"||aminoList[i][j]=="R"){
	std::pair<int,int> pr;
	pr.first = i;pr.second = j;
	posChargeList.push_back(pr);
      }
    }
  }
}


void ktlMolecule::getNegativeResidues(){
  for(int i=0;i<aminoList.size();i++){
    for(int j=0;j<aminoList[i].size();j++){
      if(aminoList[i][j]=="D"||aminoList[i][j]=="E"){
	std::pair<int,int> pr;
	pr.first = i;pr.second = j;
	negChargeList.push_back(pr);
      }
    }
  }
}


/*double ktlMolecule::getGlobalRadiusOfCurvature(){
  double globalRad =0.0;
  for(int i=0;i<hydroPhobicList.size();i++){
    double pointInf = 10000.0;
    for(int j=i+1;j<coords.size()-1;j++){
      for(int k=j+1;k<coords.size();k++){
	point p1=coords[hydroPhobicList[i].first][hydroPhobicList[i].second];
	point p2=coords[hydroPhobicList[j].first][hydroPhobicList[j].second];
	point p3=coords[hydroPhobicList[k].first][hydroPhobicList[k].second];
	double d1 = p1.eDist(p2);
	double d2 = p1.eDist(p3);
	double d3 = p2.eDist(p3);
	point dif1 =p1-p3;point dif2 = p2-p3;
	double angle = std::acos(dif1.dotprod(dif2)/(d2*d3));
	double area=2.0*std::abs(angle);
	double sphereRad = d1/area;
	if(sphereRad<pointInf){
	  pointInf=sphereRad;
	}
      }
    }
    if(i==0){
      globalRad = pointInf;
    }else{
      if(pointInf>globalRad){
	globalRad=pointInf;
      }
    }
  }
  return globalRad;
  }*/

double ktlMolecule::getGlobalRadiusOfCurvature(){
  double globalMinAv=0.0;
  for(int i=0;i<hydroPhobicList.size();i++){
    //std::cout<<i<<" "<<hydroPhobicList.size()<<" "<<hydroPhobicList[i].first<<" "<<hydroPhobicList[i].second<<"\n";
    double currMinRad = 10000.0;
    double currMinRad2 = 10000.0;
    point upperpoint,lowerpoint;
    point hydrophobPt =coords[hydroPhobicList[i].first][hydroPhobicList[i].second];
    // check if wer have the first point then need to check non neightbours
    if(hydroPhobicList[i].first==0 && hydroPhobicList[i].second==0){
      upperpoint = coords[0][1];
      point direction2 = hydrophobPt-upperpoint;direction2.normalise();
      for(int j=0;j<coords.size()-1;j++){
	for(int l=j;l<coords.size();l++){
	  for(int k=0;k<coords[j].size();k++){
	    for(int m=0;m<coords[l].size();m++){
	      point newPt =coords[j][k];
	      lowerpoint = coords[l][k];
	      point direction1 = hydrophobPt-lowerpoint;direction1.normalise();
	      // check if the
	      point direction1a = newPt-lowerpoint;direction1a.normalise();
	      point direction2a = newPt-upperpoint;direction2a.normalise();
	      double ang1 = std::acos(direction1a.dotprod(direction1));
	      double ang2 = std::acos(direction2a.dotprod(direction2));
	      if((ang1<1.57079632679) && (ang2<1.57079632679)){
		// point is inside sphere
		double d1 = lowerpoint.eDist(upperpoint);
		double d2 = lowerpoint.eDist(newPt);
		double d3 = upperpoint.eDist(newPt);
		point dif1 =lowerpoint-newPt;point dif2 = upperpoint-newPt;
		double angle = std::acos(dif1.dotprod(dif2)/(d2*d3));
		double area=2.0*std::abs(angle);
		double sphereRad = d1/area;
		if(sphereRad<currMinRad2){
		  if(sphereRad<currMinRad){
		    currMinRad = sphereRad;
		  }else{
		    currMinRad2 = sphereRad;
		  }
		}
	      }
	    }
	  }
	}
	globalMinAv = globalMinAv+currMinRad2;
      }
    }else{
      
      if(hydroPhobicList[i].second<(coords[hydroPhobicList[i].first].size()-1)){
	//upper neighbour is in the same section
	upperpoint = coords[hydroPhobicList[i].first][hydroPhobicList[i].second+1];
      }else{
	//upper neighbour first point of next curve
	upperpoint = coords[hydroPhobicList[i].first+1][0];
      }
      if(hydroPhobicList[i].second>0){
	  //lower neighbour is in the same section
	lowerpoint = coords[hydroPhobicList[i].first][hydroPhobicList[i].second-1];
      }else{
	//upper neighbour first point of next curve
	int len = coords[hydroPhobicList[i].first-1].size()-1;
	lowerpoint = coords[hydroPhobicList[i].first-1][len];
      }
      point direction1 = hydrophobPt-lowerpoint;direction1.normalise();
      point direction2 = hydrophobPt-upperpoint;direction2.normalise();
      for(int j=0;j<coords.size();j++){
	for(int k=0;k<coords[j].size();k++){
	  point newPt =coords[j][k];
	  // check if the
	  point direction1a = newPt-lowerpoint;direction1a.normalise();
	  point direction2a = newPt-upperpoint;direction2a.normalise();
	  double ang1 = std::acos(direction1a.dotprod(direction1));
	  double ang2 = std::acos(direction2a.dotprod(direction2));
	  if((ang1<1.57079632679) && (ang2<1.57079632679)){
	    // point is inside sphere
	    double d1 = lowerpoint.eDist(upperpoint);
	    double d2 = lowerpoint.eDist(newPt);
	    double d3 = upperpoint.eDist(newPt);
	    point dif1 =lowerpoint-newPt;point dif2 = upperpoint-newPt;
	    double angle = std::acos(dif1.dotprod(dif2)/(d2*d3));
	    double area=2.0*std::abs(angle);
	    double sphereRad = d1/area;
	    if(sphereRad<currMinRad2){
	      if(sphereRad<currMinRad){
		currMinRad = sphereRad;
	      }else{
		currMinRad2 = sphereRad;
	      }
	    }
	  }
	}
      }
      globalMinAv = globalMinAv+currMinRad2;
    }
  }
  return globalMinAv/hydroPhobicList.size();
}



// double ktlMolecule::getGlobalRadiusOfCurvatureBetweenSec(int &secNo){
//   double globalMinAv=0.0;
//   int firstIndex =chainList[secNo].first;
//   int secondIndex = chainList[secNo].second;
//   // get the sublist of hydrophobic residues in this section 
//   std::vector<std::pair<int,int> > hydroPhobicSubList;
//   for(int i=0;i<hydroPhobicList.size();i++){
//       if((hydroPhobicList[i].first>=firstIndex)&&(hydroPhobicList[i].first<=secondIndex)){
//         hydroPhobicSubList.push_back(hydroPhobicList[i]);
//       }
//   }
//   for(int i=0;i<hydroPhobicSubList.size();i++){
//     //std::cout<<i<<" "<<hydroPhobicList.size()<<" "<<hydroPhobicList[i].first<<" "<<hydroPhobicList[i].second<<"\n";
//     double currMinRad = 10000.0;
//     double currMinRad2 = 10000.0;
//     point upperpoint,lowerpoint;
//     point hydrophobPt =coords[hydroPhobicSubList[i].first][hydroPhobicSubList[i].second];
//     // check if we have the first or last point then need to check non neightbours
//     if((hydroPhobicSubList[i].first==0 && hydroPhobicSubList[i].second==0)||(hydroPhobicSubList[i].first==secondIndex) && (hydroPhobicSubList[i].second==(coords[secondIndex].size()-1))){
//       if(hydroPhobicSubList[i].second=(coords[hydroPhobicSubList[i].first].size()-1)){
// 	// warning this basically assumes the list has at least 2 points
// 	upperpoint = coords[hydroPhobicSubList[i].first][hydroPhobicSubList[i].second-1];
//       }else{
// 	upperpoint = coords[hydroPhobicSubList[i].first][hydroPhobicSubList[i].second+1];
//       }
//       point direction2 = hydrophobPt-upperpoint;direction2.normalise();
//       for(int j=firstIndex;j<secondIndex;j++){
// 	for(int l=j;l<=secondIndex;l++){
// 	  //check the other points are not in this section 
// 	  if((j<firstIndex && l<firstIndex)||(j<firstIndex && l>secondIndex)||(j>secondIndex && l<firstIndex)||(j>secondIndex && l>secondIndex)){
// 	    for(int k=0;k<coords[j].size();k++){
// 	      for(int m=0;m<coords[l].size();m++){
// 		point newPt =coords[j][k];
// 		lowerpoint = coords[l][k];
// 		point direction1 = hydrophobPt-lowerpoint;direction1.normalise();
// 		// check if the
// 		point direction1a = newPt-lowerpoint;direction1a.normalise();
// 		point direction2a = newPt-upperpoint;direction2a.normalise();
// 		double ang1 = std::acos(direction1a.dotprod(direction1));
// 		double ang2 = std::acos(direction2a.dotprod(direction2));
// 		if((ang1<1.57079632679) && (ang2<1.57079632679)){
// 		  // point is inside sphere
// 		  double d1 = lowerpoint.eDist(upperpoint);
// 		  double d2 = lowerpoint.eDist(newPt);
// 		  double d3 = upperpoint.eDist(newPt);
// 		  point dif1 =lowerpoint-newPt;point dif2 = upperpoint-newPt;
// 		  double angle = std::acos(dif1.dotprod(dif2)/(d2*d3));
// 		  double area=2.0*std::abs(angle);
// 		  double sphereRad = d1/area;
// 		  if(sphereRad<currMinRad2){
// 		    if(sphereRad<currMinRad){
// 		      currMinRad = sphereRad;
// 		    }else{
// 		      currMinRad2 = sphereRad;
// 		    }
// 		  }
// 		}
// 	      }
// 	    }
// 	  }
// 	}
// 	globalMinAv = globalMinAv+currMinRad2;
//       }
//     }else{
//       if(hydroPhobicSubList[i].second<(coords[hydroPhobicSubList[i].first].size()-1)){
// 	//upper neighbour is in the same section
// 	upperpoint = coords[hydroPhobicSubList[i].first][hydroPhobicSubList[i].second+1];
//       }else{
// 	//upper neighbour first point of next curve
// 	upperpoint = coords[hydroPhobicSubList[i].first+1][0];
//       }
//       if(hydroPhobicList[i].second>0){
// 	  //lower neighbour is in the same section
// 	lowerpoint = coords[hydroPhobicSubList[i].first][hydroPhobicSubList[i].second-1];
//       }else{
// 	//upper neighbour first point of next curve
// 	int len = coords[hydroPhobicSubList[i].first-1].size()-1;
// 	lowerpoint = coords[hydroPhobicSubList[i].first-1][len];
//       }
//       point direction1 = hydrophobPt-lowerpoint;direction1.normalise();
//       point direction2 = hydrophobPt-upperpoint;direction2.normalise();
//       for(int j=0;j<coords.size();j++){
// 	if((j<firstIndex)||(j>secondIndex)){
// 	  for(int k=0;k<coords[j].size();k++){
// 	    point newPt =coords[j][k];
// 	    // check if the
// 	    point direction1a = newPt-lowerpoint;direction1a.normalise();
// 	    point direction2a = newPt-upperpoint;direction2a.normalise();
// 	    double ang1 = std::acos(direction1a.dotprod(direction1));
// 	    double ang2 = std::acos(direction2a.dotprod(direction2));
// 	    if((ang1<1.57079632679) && (ang2<1.57079632679)){
// 	      // point is inside sphere
// 	      double d1 = lowerpoint.eDist(upperpoint);
// 	      double d2 = lowerpoint.eDist(newPt);
// 	      double d3 = upperpoint.eDist(newPt);
// 	      point dif1 =lowerpoint-newPt;point dif2 = upperpoint-newPt;
// 	      double angle = std::acos(dif1.dotprod(dif2)/(d2*d3));
// 	      double area=2.0*std::abs(angle);
// 	      double sphereRad = d1/area;
// 	      if(sphereRad<currMinRad2){
// 		if(sphereRad<currMinRad){
// 		  currMinRad = sphereRad;
// 		}else{
// 		  currMinRad2 = sphereRad;
// 		}
// 	      }
// 	    }
// 	  }
// 	}
//       }
//       std::cout<<" hydrophobic section "<<i<<" "<<currMinRad2<<"\n";
//       globalMinAv = globalMinAv+currMinRad2;
//     }
//   }
//   if(hydroPhobicSubList.size()>0){
//     return globalMinAv/hydroPhobicSubList.size();
//   }else{
//     return 0.0;
//   }
// }


double ktlMolecule::getGlobalRadiusOfCurvatureWithinSec(int &secNoIn,int &NoNeighbour){
  double globalMinAv=0.0;
  int secNo=secNoIn-1;
  int firstIndex =chainList[secNo].first;
  int secondIndex = chainList[secNo].second;
  // get the sublist of hydrophobic residues in this section 
  std::vector<std::pair<int,int> > hydroPhobicSubList;
  for(int i=0;i<hydroPhobicList.size();i++){
      if((hydroPhobicList[i].first>=firstIndex)&&(hydroPhobicList[i].first<=secondIndex)){
        hydroPhobicSubList.push_back(hydroPhobicList[i]);
      }
  }
  for(int i=0;i<hydroPhobicSubList.size();i++){
    std::vector<double> nearestn(NoNeighbour,10000);
    int numberNear;
    for(int j=0;j<coords.size();j++){
      point hydrophobPt =coords[hydroPhobicSubList[i].first][hydroPhobicSubList[i].second];
      if((j>=firstIndex)&& (j<= secondIndex)&&(j!=hydroPhobicSubList[i].first)){
	for(int k=0;k<coords[j].size();k++){
	  point newPt =coords[j][k];
	  double dist = newPt.eDist(hydrophobPt);
	  if(dist<nearestn[NoNeighbour-1]){
	    nearestn[NoNeighbour-1]=dist;
	    std::sort(nearestn.begin(),nearestn.end());
	  }
	}
      }
    }
    for(int i=0;i<nearestn.size();i++){
      //std::cout<<i<<" "<<nearestn[i]<<"\n";
      globalMinAv = globalMinAv+nearestn[i];
    }
  }
  if(hydroPhobicSubList.size()>0){
    return globalMinAv/double(hydroPhobicSubList.size())/double(NoNeighbour);
  }else{
    return 0.0;
  }
}


double ktlMolecule::getGlobalRadiusOfCurvatureBetweenSec(int &NoNeighbour){
  double globalRadBetween=0.0;
  int noChains = chainList.size();
  //std::cout<<"no chains are ? "<<noChains<<"\n";
  for(int i=0;i<noChains;i++){
    globalRadBetween = globalRadBetween + getGlobalRadiusOfCurvatureWithinSec(i,NoNeighbour);
  }
  return globalRadBetween;
}



    
point ktlMolecule::getCentreOfMass(std::vector<std::vector<point> > &cdSet){
  point mean(0.0,0.0,0.0);
  int noMols=0;
  for(int i=0;i<cdSet.size();i++){
    for(int j=0;j<cdSet[i].size();j++){
    mean = mean + cdSet[i][j];
    noMols=noMols+1;
    }
  }
  return mean*(1.0/double(noMols));
}

/*void ktlMolecule::rotateAndTranslate(std::vector<point> &cdSet){

  }*/


std::vector<int> ktlMolecule::checkOverlap(std::vector<std::vector<point> > &cdsIN){
  std::vector<int> overlappedSecs;
  for(int i=0;i<cdsIN.size()-1;i++){
    bool triggered =false;
    for(int j=i+1;j<cdsIN.size();j++){
      for(int k=0;k<cdsIN[i].size();k++){
	for(int l=0;l<cdsIN[j].size();l++){
	  double dist= cdsIN[i][k].eDist(cdsIN[j][l]);
	  if(k==(cdsIN[i].size()-1)&&l==0 && j==i+1){
	    dist=10.0;
	  }
	  if(dist<5.0){
	    triggered = true;
	  }
	}
      }
    }
    if(triggered==true){
      overlappedSecs.push_back(i);
    }
  }
  return overlappedSecs;
}


std::vector<double> ktlMolecule::checkOverlapWithRad(double &wRad,int &sec){
  std::vector<double> overlappedSecs;
  distSetsSecs[sec-1].clear();
  std::vector<std::vector<point> >::const_iterator firstc =coords.begin()+chainList[sec-1].first;
  std::vector<std::vector<point> >::const_iterator secondc =coords.begin()+chainList[sec-1].second+1;
  std::vector<std::vector<point> > subcoords(firstc,secondc);
  for(int i=0;i< subcoords.size()-1;i++){
    bool triggered =false;
    for(int j=i+1;j<subcoords.size();j++){
      for(int k=0;k<subcoords[i].size();k++){
	for(int l=0;l<subcoords[j].size();l++){
	  double dist= subcoords[i][k].eDist(subcoords[j][l]);
          distSetsSecs[sec-1].push_back(dist);
          if(k==(subcoords[i].size()-1)&&l==0 && j==i+1){
	    dist=10.0;
          }
	  if(dist<wRad){
	    triggered = true;
	    overlappedSecs.push_back(dist);
	  }
	}
      }
    }
    //if(triggered==true){
    // overlappedSecs.push_back(i);
    //}
  }
  return overlappedSecs;
}

std::vector<double> ktlMolecule::checkOverlapWithRad(double &wRad){
  std::vector<double> overlappedSecs;
  distSets.clear();
  for(int i=0;i<coords.size();i++){
    bool triggered =false;
    for(int j=i;j<coords.size();j++){
      for(int k=0;k<coords[i].size();k++){
	for(int l=0;l<coords[j].size();l++){
	  double dist= coords[i][k].eDist(coords[j][l]);
	  if(i==j){
	    if(l>k){
	      distSets.push_back(dist);
	    }
	  }else{
	    distSets.push_back(dist);
	  }
          if((k==(coords[i].size()-1)&&l==0 && j==i+1)||j==i){
	    dist=10.0;
          }
	  if(dist<wRad){
	    triggered = true;
	    overlappedSecs.push_back(dist);
	  }
	}
      }
    }
    //if(triggered==true){
    // overlappedSecs.push_back(i);
    //}
  }
  return overlappedSecs;
}


std::vector<double> ktlMolecule::getDistSet(){
  return distSets;
};

double ktlMolecule::compareDistances(std::vector<std::vector<point> > &coords2){
  std::vector<int> overlappedSecs;
  double distDifTot =0.0;
  double n=0;;
  for(int i=0;i<coords.size()-1;i++){
    for(int j=i+1;j<coords.size();j++){
      for(int k=0;k<coords[i].size();k++){
	for(int l=0;l<coords[j].size();l++){
	  double dist= coords[i][k].eDist(coords[j][l]);
	  double dist2= coords2[i][k].eDist(coords2[j][l]);
	  distDifTot = distDifTot+ std::abs(dist2-dist);
	  n=n+1.0;
	}
      }
    }
  }
  return distDifTot/n;
}


point ktlMolecule::binorm(point &T1,point &T2){
  point cp = T1.cross(T2);
  double cpsize = cp.length();
  if(cpsize>0.0000000001){
    cp.normalise();
  }else{
    point vert(0.0,0.0,1.0);
    if(std::abs(T1.dotprod(vert))<0.99999999){
      point vert(0.0,0.0,1.0);
       point cp = vert.cross(T1);
       cp.normalise();
    }else{
      point vert(0.0,1.0,0.0);
      cp = vert.cross(T1);
      cp.normalise();
    }
  }
  return cp;
}

point ktlMolecule::parallelTransport(point &tan1,point &tan2,point &norm1){
  point bvec = binorm(tan1,tan2);
  bvec.normalise();
  if(std::abs(tan1.dotprod(tan2))<0.9999999){
    double costhe = tan1.dotprod(tan2);
    double sinthe = std::sqrt(1-costhe*costhe);
    return norm1*costhe + (bvec.cross(norm1))*sinthe + bvec*(bvec.dotprod(norm1))*(1.0-costhe);
  }else{
    return norm1;
  }
}

std::vector<std::vector<point> > ktlMolecule::updateFrame(std::vector<point> &section,point &tangent,point &normal,point &binormal){
  point tanOld = tangent;
  point currNorm = normal;
  point newNorm,newBinorm;
  std::vector<std::vector<point> > frameSet;
  for(int i=1;i<section.size();i++){
    point tanNew = section[i]-section[i-1];
    tanNew.normalise();
    newNorm = parallelTransport(tanOld,tanNew,currNorm);
    currNorm = newNorm;
    newBinorm = tanNew.cross(currNorm);
    std::vector<point> frame;
    frame.push_back(tanNew);frame.push_back(currNorm);frame.push_back(newBinorm);
    frameSet.push_back(frame);
    tanOld = tanNew;
  }
  return frameSet;
}



void ktlMolecule::getFrameForBackbone(){
  for(int i=0;i<chainList.size();i++){
    //std::cout<<"in loop ? "<<i<<" "<<chainList[i].first<<" "<<chainList[i].second<<"\n";
    std::vector<std::vector<point> >::const_iterator first;
    std::vector<std::vector<point> >::const_iterator second;
    if(i==0){
      first =coords.begin()+chainList[i].first;
      second =coords.begin()+chainList[i].second;
    }else{
      first =coords.begin()+chainList[i].first;
      second =coords.begin()+chainList[i].second;
    }
    std::vector<std::vector<point> > subsec(first,second);
    point currTan = subsec[0][1]-subsec[0][0];
    currTan.normalise();
    point vert(0.0,0.0,1.0);
    point currNor;
    if(std::abs(currTan.dotprod(vert))<0.99999999){
      currNor = vert.cross(currTan);
      currNor.normalise();
    }else{
      point vert2(0.0,1.0,0.0);
      currNor = vert2.cross(currTan);
      currNor.normalise();
    }
    point currBi = currTan.cross(currNor);
    for(int i=0;i<subsec.size();i++){
      //std::cout<<i<<" "<<" "<<subsec[i].size()<<" "<<subsec.size()<<"\n";
      std::vector<std::vector<point> > newFrms =updateFrame(subsec[i],currTan,currNor,currBi);
      int sz = newFrms.size()-1;
      std::vector<point> tnsec;
      std::vector<point> nmsec;
      std::vector<point> bnsec;
      for(int j=0;j<=sz;j++){
	tnsec.push_back(newFrms[j][0]);
	nmsec.push_back(newFrms[j][1]);
	bnsec.push_back(newFrms[j][2]);
      }
      tanlist.push_back(tnsec);
      normlist.push_back(nmsec);
      binormlist.push_back(bnsec);
      currTan = newFrms[sz][0];currNor = newFrms[sz][1];currBi = newFrms[sz][2];
    }
  }
}

bool ktlMolecule::checkCalphas(std::vector<std::vector<point> > &coordsIn){
  bool violated=false;
  for(int i=0;i<coordsIn.size();i++){
    for(int j=0;j<coordsIn[i].size();j++){
      double dist =0.0;
      if((j==coordsIn[i].size()-1) && (i<(coordsIn.size()-1))){
	dist = coordsIn[i][j].eDist(coordsIn[i+1][0]);
      }else if((j==coordsIn[i].size()-1) && (i==coordsIn.size()-1)){
	dist = 3.7;
      }
      else{
	dist = coordsIn[i][j].eDist(coordsIn[i][j+1]);
      }
      if((dist>4.0)|| (dist<3.3)){
	violated =true;
      }
    }
  }
  return violated;
}

bool ktlMolecule::checkCalphas(){
  bool violated=false;
  for(int i=0;i<coords.size();i++){
    for(int j=0;j<coords[i].size();j++){
      double dist =0.0;
      if((j==coords[i].size()-1) && (i<(coords.size()-1))){
	dist = coords[i][j].eDist(coords[i+1][0]);
      }else if((j==coords[i].size()-1) && (i==coords.size()-1)){
	dist = 3.7;
      }
      else{
	dist = coords[i][j].eDist(coords[i][j+1]);
      }
      if((dist>4.0)|| (dist<3.3)){
	violated =true;
      }
    }
  }
  return violated;
}


bool ktlMolecule::checkCalphas(int &secNo){
  int sec = secNo-1;
  std::vector<std::vector<point> >::const_iterator firstc =coords.begin()+chainList[sec].first;
  std::vector<std::vector<point> >::const_iterator secondc =coords.begin()+chainList[sec].second+1;
  std::vector<std::vector<point> > coordsSec(firstc,secondc);
  bool violated=false;
  for(int i=0;i<coordsSec.size();i++){
    //std::cout<<coordsSec[i].size()<<"\n";
    for(int j=0;j<coordsSec[i].size();j++){
      double dist =0.0;
      if((j==coordsSec[i].size()-1) && (i<(coordsSec.size()-1))){
	dist = coordsSec[i][j].eDist(coordsSec[i+1][0]);
      }else if((j==coordsSec[i].size()-1) && (i==coordsSec.size()-1)){
	dist = 3.7;
      }
      else{
	dist = coordsSec[i][j].eDist(coordsSec[i][j+1]);
      }
      if((dist>4.3)|| (dist<2.6)){
	//std::cout<<i<<" "<<j<<" "<<dist<<"\n";
	violated =true;
      }
    }
  }
  return violated;
}



int ktlMolecule::getRandomMolecule(){
  coords.clear();
  int noOverLaps=0;
  for(int i=0;i<chainList.size();i++){
     std::cout<<"chain "<<i<<"\n";
    std::vector<std::pair<std::string,int> >::const_iterator first =nameSizeList.begin()+chainList[i].first;
    std::vector<std::pair<std::string,int> >::const_iterator second =nameSizeList.begin()+chainList[i].second+1;
    std::vector<std::pair<std::string,int> > subns(first,second);
    point sp(0.0,0.0,0.0);
    bool suc=true;
    std::vector<std::vector<point> > submol = rmg.makeRandomMolecule(subns,sp,suc);
    std::cout<<"made man\n?";
    std::vector<int> overLapList = checkOverlap(submol);
    std::vector<int> currOverLapList= overLapList;
    std::random_device rdev{};
    std::default_random_engine generator1{rdev()};
    int l=0;
    while(l<10000 && currOverLapList.size()>0){
      std::vector<std::vector<point> > tempCoords = submol;
      //std::cout<<currOverLapList.size()<<" "<<tempCoords.size()<<" "<<l<<"\n";
      if(currOverLapList.size()>l){
	std::uniform_int_distribution<int> distributionR(0,currOverLapList.size()-1);
	int iv= distributionR(generator1);
	if(currOverLapList[iv]==0){
	  changeMoleculeSingle(currOverLapList[1],tempCoords,subns);  
	}else{
	  changeMoleculeSingle(currOverLapList[iv],tempCoords,subns);  
	}
      }else{
	changeMoleculeSingle(currOverLapList[0],tempCoords,subns);
      }
      overLapList = checkOverlap(tempCoords);
      bool caca = checkCalphas(tempCoords);
      if(overLapList.size()<currOverLapList.size() && caca==false){
	currOverLapList=overLapList;
	submol = tempCoords;
      }
      l++;
    }
    noOverLaps = noOverLaps + currOverLapList.size();
    // now add the section to the whole
    for(int j=0;j<submol.size();j++){
      coords.push_back(submol[j]);
    } 
  }
  return noOverLaps;
}

/*int ktlMolecule::getRandomMoleculeReset(){
   coords.clear();
  for(int i=0;i<chainList.size();i++){
    std::vector<std::pair<std::string,int> >::const_iterator first =nameSizeList.begin()+chainList[i].first;
    std::vector<std::pair<std::string,int> >::const_iterator second =nameSizeList.begin()+chainList[i].second+1;
    std::vector<std::pair<std::string,int> > subns(first,second);
    point sp(0.0,0.0,0.0);
    bool suc=true;
    std::vector<std::vector<point> > submol = rmg.makeRandomMolecule(subns,sp,suc);
    for(int j=0;j<submol.size();j++){
      coords[j]=submol[j];
    } 
  }
  std::vector<int> overLapList = checkOverlap(coords);
  std::vector<int> currOverLapList= overLapList;
  std::random_device rdev{};
  std::default_random_engine generator1{rdev()};
  int l=0;
  while(l<10000 && currOverLapList.size()>0){
    std::vector<std::vector<point> > tempCoords = coords;
    //std::cout<<currOverLapList.size()<<"\n";
    if(currOverLapList.size()>l){
      std::uniform_int_distribution<int> distributionR(0,currOverLapList.size()-1);
      int iv= distributionR(generator1);
      if(currOverLapList[iv]==0){
	changeMoleculeSingle(currOverLapList[1],tempCoords);  
      }else{
	changeMoleculeSingle(currOverLapList[iv],tempCoords);  
      }
    }else{
      changeMoleculeSingle(currOverLapList[0],tempCoords);
    }
    overLapList = checkOverlap(tempCoords);
    if(overLapList.size()<currOverLapList.size()){
      currOverLapList=overLapList;
      coords = tempCoords;
    }
    l++;
  }
  return currOverLapList.size();
  }*/


void ktlMolecule::getRandomMoleculeAllowOverlap(){
  for(int i=0;i<chainList.size();i++){
    std::vector<std::pair<std::string,int> >::const_iterator first =nameSizeList.begin()+chainList[i].first;
    std::vector<std::pair<std::string,int> >::const_iterator second =nameSizeList.begin()+chainList[i].second+1;
    std::vector<std::pair<std::string,int> > subns(first,second);
    point sp(0.0,0.0,0.0);
    bool suc=true;
    std::vector<std::vector<point> > submol = rmg.makeRandomMolecule(subns,sp,suc);
    for(int j=0;j<submol.size();j++){
      coords.push_back(submol[j]);
    } 
  }
}


void ktlMolecule::getRandomMoleculeAllowOverlapReset(){
  for(int i=0;i<chainList.size();i++){
    std::vector<std::pair<std::string,int> >::const_iterator first =nameSizeList.begin()+chainList[i].first;
    std::vector<std::pair<std::string,int> >::const_iterator second =nameSizeList.begin()+chainList[i].second+1;
    std::vector<std::pair<std::string,int> > subns(first,second);
    point sp(0.0,0.0,0.0);
    bool suc=true;
    std::vector<std::vector<point> > submol = rmg.makeRandomMolecule(subns,sp,suc);
    for(int j=0;j<submol.size();j++){
      coords[j]=submol[j];
    } 
  }
}

void ktlMolecule::removeOverlap(){
  std::vector<int> overLapList = checkOverlap(coords);
  std::vector<int> currOverLapList= overLapList;
  std::random_device rdev{};
  std::default_random_engine generator1{rdev()};
  int l=0;
  std::vector<std::pair<std::string,int> > nameSizeListTemp = nameSizeList;
  while(l<10000 && currOverLapList.size()>0){
    std::vector<std::vector<point> > tempCoords = coords;
    //std::cout<<currOverLapList.size()<<"\n";
    if(currOverLapList.size()>l){
      std::uniform_int_distribution<int> distributionR(0,currOverLapList.size()-1);
      int iv= distributionR(generator1);
      if(currOverLapList[iv]==0){
	changeMoleculeSingle(currOverLapList[1],tempCoords,nameSizeListTemp);  
      }else{
	changeMoleculeSingle(currOverLapList[iv],tempCoords,nameSizeListTemp);  
      }
    }else{
      changeMoleculeSingle(currOverLapList[0],tempCoords,nameSizeListTemp);
    }
    overLapList = checkOverlap(tempCoords);
    if(overLapList.size()<currOverLapList.size()){
      currOverLapList=overLapList;
      coords = tempCoords;
    }
    l++;
  }
}


/*void ktlMolecule::resetRandomMolecule(){
  for(int i=0;i<chainList.size();i++){
    std::vector<std::pair<std::string,int> >::const_iterator first =nameSizeList.begin()+chainList[i].first;
    std::vector<std::pair<std::string,int> >::const_iterator second =nameSizeList.begin()+chainList[i].second+1;
    std::vector<std::pair<std::string,int> > subns(first,second);
    point sp(0.0,0.0,0.0);
    bool suc=true;
    std::vector<std::vector<point> > submol = rmg.makeRandomMolecule(subns,sp,suc);
    for(int j=0;j<submol.size();j++){
      coords[chainList[i].second+1+j]=submol[j];
    } 
  }
  std::vector<int> overLapList = checkOverlap(coords);
  std::vector<int> currOverLapList= overLapList;
  std::random_device rdev{};
  std::default_random_engine generator1{rdev()};
  int l=0;
  while(l<10000 && currOverLapList.size()>0){
    std::vector<std::vector<point> > tempCoords = coords;
    if(currOverLapList.size()>1){
      std::uniform_int_distribution<int> distributionR(0,currOverLapList.size()-1);
      int iv= distributionR(generator1);
      if(currOverLapList[iv]==0){
	changeMoleculeSingle(currOverLapList[1],tempCoords,);  
      }else{
	changeMoleculeSingle(currOverLapList[iv],tempCoords);  
      }
    }else{
      changeMoleculeSingle(currOverLapList[0],tempCoords);
    }
    overLapList = checkOverlap(tempCoords);
    if(overLapList.size()<currOverLapList.size()){
      currOverLapList=overLapList;
      coords = tempCoords;
    }
    l++;
    std::cout<<l<<" "<<overLapList.size()<<"\n";
  }
  }*/

void ktlMolecule::changeMoleculeSingle(int &index){
  point sp(0.0,0.0,0.0);
  //std::cout<<index<<" "<<coords[index].size()<<"\n";
  bool suc=true;
  if(index<coords.size()-2){
    // check for repeat helices
    if(nameSizeList[index].first=="loop"&&nameSizeList[index+1].first=="Helix"&&nameSizeList[index+2].first=="Helix"){
      //change loo then encounter repeat helices
       //change helix and encounter repreat helices
      int noHelices=0;int k=index+1;
      while(k<coords.size()&&nameSizeList[k].first=="Helix"){
	noHelices++;
	k++;
      }
      maxDistChange =rmg.reshapeMolLoopThenHelixSet(nameSizeList,coords,index,noHelices,sp,suc);
    }else if(nameSizeList[index].first=="Helix"&&nameSizeList[index+1].first=="Helix"&&nameSizeList[index+2].first=="Helix"){
      //change helix and encounter repreat helices
      int noHelices=0;int k=index+1;
      while(k<coords.size()&&nameSizeList[k].first=="Helix"){
	noHelices++;
	k++;
      }
      maxDistChange =rmg.reshapeMolHelixSet(nameSizeList,coords,index,noHelices,sp,suc);
    }else{
      maxDistChange =rmg.reshapeMol(nameSizeList,coords,index,sp,suc);
    }
  }else{
    maxDistChange =rmg.reshapeMol(nameSizeList,coords,index,sp,suc);
  }
}

void ktlMolecule::changeMoleculeSet(std::vector<int> &indicies){
  point sp(0.0,0.0,0.0);
  //std::cout<<index<<" "<<coords[index].size()<<"\n";
  bool suc=true;
  for(int i=0;i<indicies.size();i++){
    int index = indicies[i];
    if(index<coords.size()-2){
    // check for repeat helices
    if(nameSizeList[index].first=="loop"&&nameSizeList[index+1].first=="Helix"&&nameSizeList[index+2].first=="Helix"){
      //change loo then encounter repeat helices
       //change helix and encounter repreat helices
      int noHelices=0;int k=index+1;
      while(k<coords.size()&&nameSizeList[k].first=="Helix"){
	noHelices++;
	k++;
      }
      maxDistChange =rmg.reshapeMolLoopThenHelixSet(nameSizeList,coords,index,noHelices,sp,suc);
    }else if(nameSizeList[index].first=="Helix"&&nameSizeList[index+1].first=="Helix"&&nameSizeList[index+2].first=="Helix"){
      //change helix and encounter repreat helices
      int noHelices=0;int k=index+1;
      while(k<coords.size()&&nameSizeList[k].first=="Helix"){
	noHelices++;
	k++;
      }
      maxDistChange =rmg.reshapeMolHelixSet(nameSizeList,coords,index,noHelices,sp,suc);
    }else{
      maxDistChange =rmg.reshapeMol(nameSizeList,coords,index,sp,suc);
    }
  }else{
    maxDistChange =rmg.reshapeMol(nameSizeList,coords,index,sp,suc);
  }
  }
}


void ktlMolecule::changeMoleculeSingleMulti(int &index,int secIn){
  int sec = secIn-1;
  point sp(0.0,0.0,0.0);
  std::vector<std::pair<std::string,int> >::const_iterator first =nameSizeList.begin()+chainList[sec].first;
  std::vector<std::pair<std::string,int> >::const_iterator second =nameSizeList.begin()+chainList[sec].second+1;
  std::vector<std::pair<std::string,int> > subns(first,second);
  std::vector<std::vector<point> >::const_iterator firstc =coords.begin()+chainList[sec].first;
  std::vector<std::vector<point> >::const_iterator secondc =coords.begin()+chainList[sec].second+1;
  std::vector<std::vector<point> > subcoords(firstc,secondc);
  bool suc=true;
  if(index<subcoords.size()-2){
    // // check for repeat heliced
    if(subns[index].first=="loop"&&subns[index+1].first=="Helix"&&subns[index+2].first=="Helix"){
      //change loo then encounter repeat helices
       //change helix and encounter repreat helices
      int noHelices=0;int k=index+1;
      while(k<subcoords.size()&&subns[k].first=="Helix"){
	noHelices++;
	k++;
      }
      maxDistChange =rmg.reshapeMolLoopThenHelixSet(subns,subcoords,index,noHelices,sp,suc);
    }else if(subns[index].first=="Helix"&&subns[index+1].first=="Helix"&&subns[index+2].first=="Helix"){
      //change helix and encounter repreat helices
      int noHelices=0;int k=index+1;
      while(k<subcoords.size()&&subns[k].first=="Helix"){
	noHelices++;
	k++;
      }
      maxDistChange =rmg.reshapeMolHelixSet(subns,subcoords,index,noHelices,sp,suc);
    }else{
      maxDistChange =rmg.reshapeMol(subns,subcoords,index,sp,suc);
    }
  }else{
    maxDistChange =rmg.reshapeMol(subns,subcoords,index,sp,suc);
  }
  // now replace the new sections into the main chain
  int d = chainList[sec].second+1-chainList[sec].first;
  std::copy_n(subcoords.begin(),d,&coords[chainList[sec].first]);
}

void ktlMolecule::changeMoleculeSetMulti(std::vector<int>  &indicies,int secIn){
  point sp(0.0,0.0,0.0);
  int sec = secIn-1;
  std::vector<std::pair<std::string,int> >::const_iterator first =nameSizeList.begin()+chainList[sec].first;
  std::vector<std::pair<std::string,int> >::const_iterator second =nameSizeList.begin()+chainList[sec].second+1;
  std::vector<std::pair<std::string,int> > subns(first,second);
  std::vector<std::vector<point> >::const_iterator firstc =coords.begin()+chainList[sec].first;
  std::vector<std::vector<point> >::const_iterator secondc =coords.begin()+chainList[sec].second+1;
  std::vector<std::vector<point> > subcoords(firstc,secondc);
  bool suc=true;
  for(int i=0;i<indicies.size();i++){
    int index=indicies[i];
    if(index<subcoords.size()-2){
      // check for repeat helices
      if(subns[index].first=="loop"&&subns[index+1].first=="Helix"&&subns[index+2].first=="Helix"){
	//change loo then encounter repeat helices
	//change helix and encounter repreat helices
	int noHelices=0;int k=index+1;
	while(k<subcoords.size()&&subns[k].first=="Helix"){
	  noHelices++;
	  k++;
	}
	maxDistChange =rmg.reshapeMolLoopThenHelixSet(subns,subcoords,index,noHelices,sp,suc);
      }else if(subns[index].first=="Helix"&&subns[index+1].first=="Helix"&&subns[index+2].first=="Helix"){
	//change helix and encounter repreat helices
	int noHelices=0;int k=index+1;
	while(k<subcoords.size()&&subns[k].first=="Helix"){
	  noHelices++;
	k++;
	}
	maxDistChange =rmg.reshapeMolHelixSet(subns,subcoords,index,noHelices,sp,suc);
      }else{
	maxDistChange =rmg.reshapeMol(subns,subcoords,index,sp,suc);
      }
    }else{
      maxDistChange =rmg.reshapeMol(subns,subcoords,index,sp,suc);
    }
  }
  // now replace the new sections into the main chain
  int d = chainList[sec].second+1-chainList[sec].first;
  std::copy_n(subcoords.begin(),d,&coords[chainList[sec].first]);
}

void ktlMolecule::changeMoleculeMultiRotate(double &angle,point &k,int secIn,point &transVec){
  point sp(0.0,0.0,0.0);
  int sec = secIn-1;
  std::vector<std::vector<point> >::const_iterator firstc =coords.begin()+chainList[sec].first;
  std::vector<std::vector<point> >::const_iterator secondc =coords.begin()+chainList[sec].second+1;
  std::vector<std::vector<point> > subcoords(firstc,secondc);
  bool suc=true;
  point com = getCentreOfMass(subcoords);
  rotateSection(subcoords,com,k,angle,transVec);
  // now replace the new sections into the main chain
  int di = chainList[sec].second+1-chainList[sec].first;
  std::copy_n(subcoords.begin(),di,&coords[chainList[sec].first]);
}

void ktlMolecule::replicateMolecule(int &noReplications){
  std::random_device rdev{};
  std::default_random_engine generator1{rdev()};
  std::uniform_real_distribution<> rotAng(0.0,0.5);
  std::uniform_real_distribution<> theAng(0.0,3.14159265359);
  std::uniform_real_distribution<> phiAng(0.0,6.28318530718); 
  point com = getCentreOfMass(coords);
  double maxRad=0.0;
  for(int i=0;i<coords.size();i++){
    for(int j=0;j<coords[i].size();j++){
      double rad = com.eDist(coords[i][j]);
      if(rad>maxRad){
        maxRad = rad;
      }
    }
  } 
  std::vector<std::vector<point> > replicatedCoords = coords;
  std::vector<std::pair<std::string,int> > replicatedNameSizeList=nameSizeList;
  std::vector<std::pair<int,int> > chainListCopy = chainList;
  std::vector<std::vector<std::string> > aminoListCopy = aminoList;
  for(int i=1;i<=noReplications;i++){
    std::vector<std::vector<point> > newPts = coords;
    double angle = rotAng(generator1);double theta=theAng(generator1);double phi = phiAng(generator1);
    point k(std::sin(theta)*std::cos(phi),std::sin(theta)*std::sin(phi),std::cos(theta));
    theta=theAng(generator1);phi = phiAng(generator1);
    point translate(std::sin(theta)*std::cos(phi),std::sin(theta)*std::sin(phi),std::cos(theta));
    translate.scalarMult(2.2*maxRad);
    rotateSection(newPts,com,k,angle,translate);
    // now append
    replicatedCoords.insert(replicatedCoords.end(),newPts.begin(),newPts.end());
    // now add on to the nameList separator 
    std::vector<std::pair<std::string,int> > copyNameList;
    for(int j=0;j<nameSizeList.size();j++){
      std::pair<std::string,int> pr;  
      pr.first=nameSizeList[j].first;
      pr.second=nameSizeList[j].second;
      copyNameList.push_back(pr);
    }
    replicatedNameSizeList.insert(replicatedNameSizeList.end(),copyNameList.begin(),copyNameList.end());
    std::pair<int,int> prInts;
    prInts.first=chainListCopy[chainListCopy.size()-1].second+1;
    prInts.second=chainListCopy[chainListCopy.size()-1].second+chainList[chainList.size()-1].second+1;
    chainListCopy.push_back(prInts);
    //finally replicate the aminoList and hydration list
    aminoListCopy.insert(aminoListCopy.begin(),aminoList.begin(),aminoList.end());
  }  
  aminoList = aminoListCopy;
  nameSizeList=replicatedNameSizeList;
  chainList=chainListCopy;
  coords=replicatedCoords;
}


void ktlMolecule::changeMoleculeLocal(int &index,double variationSize){
  point sp(0.0,0.0,0.0);
  //std::cout<<index<<" "<<coords[index].size()<<"\n";
  bool suc=true;
  maxDistChange = rmg.reshapeMolSmallVariation(nameSizeList,coords,index,sp,variationSize,suc);
}

void ktlMolecule::changeMoleculeLocalSet(std::vector<int> &indicies,double variationSize){
  point sp(0.0,0.0,0.0);
  //std::cout<<index<<" "<<coords[index].size()<<"\n";
  bool suc=true;
  for(int i=0;i<indicies.size();i++){
    maxDistChange = rmg.reshapeMolSmallVariation(nameSizeList,coords,indicies[i],sp,variationSize,suc);
  }
}

/*void ktlMolecule::changeMoleculeSingle(int &index,std::vector<std::vector<point> > &cdsIn,std::vector<std::pair<std::string,int> > &nameSizeSubList){
  point sp(0.0,0.0,0.0);
  bool suc=true;
  if(index<cdsIn.size()-2){
    // check for repeat helices
    if(nameSizeSubList[index].first=="loop"&&nameSizeSubList[index+1].first=="Helix"&&nameSizeSubList[index+2].first=="Helix"){
      //change loo then encounter repeat helices
       //change helix and encounter repreat helices
      int noHelices=0;int k=index+1;
      while(k<coords.size()&&nameSizeSubList[k].first=="Helix"){
	noHelices++;
	k++;
      }
      std::cout<<"here hel hel\n";
      maxDistChange =rmg.reshapeMolLoopThenHelixSet(nameSizeSubList,cdsIn,index,noHelices,sp,suc);
    }else if(nameSizeSubList[index+1].first=="Helix"&&nameSizeSubList[index+2].first=="Helix"){
      //change helix and encounter repreat helices
      int noHelices=0;int k=index+1;
      while(k<cdsIn.size()&&nameSizeSubList[k].first=="Helix"){
	noHelices++;
	k++;
      }
      std::cout<<"change hel set\n";
      maxDistChange =rmg.reshapeMolHelixSet(nameSizeSubList,cdsIn,index,noHelices,sp,suc);
    }else{
      std::cout<<"change normal 1\n";
      maxDistChange =rmg.reshapeMol(nameSizeSubList,cdsIn,index,sp,suc);
    }
  }else{
    std::cout<<"change normal 2\n";
    maxDistChange =rmg.reshapeMol(nameSizeSubList,cdsIn,index,sp,suc);
  }
}
*/

void ktlMolecule::changeMoleculeSingle(int &index,std::vector<std::vector<point> > &cdsIn,std::vector<std::pair<std::string,int> > &nameSizeSubList){
  point sp(0.0,0.0,0.0);
  bool suc=true;
  maxDistChange =rmg.reshapeMol(nameSizeSubList,cdsIn,index,sp,suc);
}

void ktlMolecule::rotation3D(point &p,point &centre,point &k,double &cosangle,double &sinangle){
  point out;
  // translate so the centre is at the origin
  out = p-centre;
  point term1 = out*cosangle;
  point term2 = k.cross(out);term2 = term2*sinangle;
  double dotProd = k.dotprod(out);
  point term3 = k*dotProd;
  term3 = term3*(1.0-cosangle);
  out = term1+term2 +term3;
  out = out +centre;
  p=out;
}

void ktlMolecule::rotateSection(std::vector<std::vector<point> >  &section,point &centre,point &k,double &angle,point &transVec){
  double cosangle = std::cos(angle);double sinangle = std::sin(angle);
  for(int i=0;i<section.size();i++){
    for(int j=0;j<section[i].size();j++){
      rotation3D(section[i][j],centre,k,cosangle,sinangle);
      section[i][j] = section[i][j]+transVec;
    }
  }
}


bool ktlMolecule::checkOverlapSbond(double &minDist){
  int mSize = coords.size();
  int i=0;int j=0;
  bool exceededSz = false;
  while(i<mSize-1 && exceededSz== false){
    j=i+1;
    while(j<mSize  && exceededSz== false){
      //looking at section i and j
      for(int k=0;k<coords[i].size();k++){
	for(int l=0;l<coords[j].size();l++){
	  double dst = coords[i][k].eDist(coords[j][l]);
	  if(dst<minDist){
	    if(j!= i+1){
	       exceededSz=true;
	    }else{
	      // make sure this is not the nearest neighbour crossing
	      if(l==0 && k==(coords[i].size()-1)){
		exceededSz=false;
	      }else{
		//std::cout<<"here "<<i<<" "<<j<<" "<<l<<" "<<k<<" "<<(coords[i].size()-1)<<" "<<dst<<"\n";
		exceededSz=true;
	      }
	    }
	  }
	}
      }
      j++;
    }
    i++;
  }
  //std::cout<<exceededSz<<"\n";
  return exceededSz;
}


void ktlMolecule::changeMoleculeSingleCheckOverlap(){
  point sp(0.0,0.0,0.0);
  bool suc=true;
  bool isDone=false;
  std::random_device rdev{};
  std::default_random_engine generator1{rdev()};
  std::uniform_int_distribution<int> distributionR(2,coords.size()-2);
  while(isDone==false){
    std::vector<std::vector<point> > tempCoords = coords;
    int iv= distributionR(generator1);
    maxDistChange =rmg.reshapeMol(nameSizeList,tempCoords,iv,sp,suc);
    std::vector<int> overLapList = checkOverlap(tempCoords);
    std::cout<<overLapList.size()<<"\n";
    if(overLapList.size()==0){
      isDone=true;
      std::cout<<iv<<" hmm ? "<<coords.size()<<"\n";
      coords = tempCoords;
    }
  }
}

void ktlMolecule::writeMoleculeToFile(const char* filename){
  std::ofstream ofile;
  ofile.open(filename);
  //std::cout<<filename<<"\n";
  if(ofile.is_open()){
     for(int k=0;k<chainList.size();k++){
      for(int i=chainList[k].first;i<=chainList[k].second;i++){
       for(int j=0;j<coords[i].size();j++){
	 //std::cout<<coords[i][j].getX()<<" "<<coords[i][j].getY()<<" "<<coords[i][j].getZ()<<"\n"; 
	   ofile<<coords[i][j].getX()<<" "<<coords[i][j].getY()<<" "<<coords[i][j].getZ()<<"\n"; 
       }
       ofile<<"\n";
      }
      ofile<<"End chain "<<k+1<<"\n";
    }
  }else{
    ofile<<"cannot write molecule to file";
  }
  ofile.close();
}

void ktlMolecule::getBackboneStats(){
  std::cout<<"how big "<<coords.size()<<"\n";
  // get the pair pair distances, check something isn't wrong 
  point currCoord=coords[0][1];point prevCoord=coords[0][0];
  int section =0;int subsection=0;
  double maxdist = 0.0;int joinOrInterior=0;
  int npts =1;
  for(int i=0;i<coords.size();i++){
    if(i==0){
      for(int j=1;j<coords[i].size();j++){
        npts++;
        currCoord = coords[i][j];
        double dist= currCoord.eDist(prevCoord);
        if(dist>maxdist){
         maxdist = dist;
         section = i;
         subsection =j;  
          if(j==0){
            joinOrInterior =1;
          }else{
            joinOrInterior =0;
          }
        }
        prevCoord=currCoord;
      }
    }
    else{
       for(int j=0;j<coords[i].size();j++){
        npts++;
        currCoord = coords[i][j];
        double dist= currCoord.eDist(prevCoord);
        if(dist>maxdist){
         maxdist = dist; 
         section = i;  
         subsection =j;  
          if(j==0){
            joinOrInterior =1;
          }else{
            joinOrInterior =0;
          }
        }
        prevCoord=currCoord;
      }
    }
  }
  if(joinOrInterior==1){
    std::cout<<"biggest distance on link "<<maxdist<<"\n";
    std::cout<<"number of points is "<<npts<<"\n";
    std::cout<<"it occured with end point at "<<section<<" point "<<subsection<<" "<<coords[section][0].eDist(coords[section-1][coords[section-1].size()-1])<<"\n";
    std::cout<<"which is a join from "<<nameSizeList[section-1].first<<" "<<nameSizeList[section].first<<"\n";
  }else{
    std::cout<<"biggest distance on interior "<<maxdist<<"\n";
    std::cout<<"number of points is "<<npts<<"\n";
    std::cout<<"it occured with end point at "<<section<<" point "<<subsection<<"\n";
    std::cout<<"which is on interior "<<nameSizeList[section].first<<"\n";
  }
}

std::vector<std::pair<double,double> > ktlMolecule::getKapTauVals(){
  std::vector<std::pair<double,double> > ktllst;
  for(int k=0;k<chainList.size();k++){
    for(int i=chainList[k].first;i<=chainList[k].second;i++){
      for(int j=0;j<coords[i].size();j++){
         if(i!=chainList[k].second){
          if(j<coords[i].size()-3){
            ktllst.push_back(rmg.kapTau(coords[i][j],coords[i][j+1],coords[i][j+2],coords[i][j+3]));
          }else if(j==coords[i].size()-3){
            ktllst.push_back(rmg.kapTau(coords[i][j],coords[i][j+1],coords[i][j+2],coords[i+1][0]));
          }else if(j==coords[i].size()-2){
            ktllst.push_back(rmg.kapTau(coords[i][j],coords[i][j+1],coords[i+1][0],coords[i+1][1]));
          }else if(j==coords[i].size()-1){
            ktllst.push_back(rmg.kapTau(coords[i][j],coords[i][0],coords[i+1][1],coords[i+1][2]));
          }
        }else{
          if(j<coords[i].size()-3){
            ktllst.push_back(rmg.kapTau(coords[i][j],coords[i][j+1],coords[i][j+2],coords[i][j+3]));
          }
        }
      }
    }
  }
  return ktllst; 
}


double ktlMolecule::getPairDistance(std::pair<int,int> &index1,std::pair<int,int> &index2){
  point pt1 = coords[index1.first][index1.second];
  point pt2 = coords[index2.first][index2.second];
  //std::cout<<index1.first<<" "<<index1.second<<"\n";
  //pt1.printPoint();
  //std::cout<<index2.first<<" "<<index2.second<<"\n";
  //pt2.printPoint();
  return pt1.eDist(pt2);
}

double ktlMolecule::lennardJones(double &idealDist,double &currDist,int noConPred,double &weightCoeff){
  //for actual Leennard jones
  /*double rat = idealDist/currDist;
  double edist = 1.0;
  double sixpow =rat*rat*rat*rat*rat*rat;
  return (edist*(sixpow*sixpow-2.0*sixpow)+edist)/double(noConPred);
  */
  double dif = (idealDist-currDist)/idealDist;
  return weightCoeff*dif*dif;
}

double ktlMolecule::getFitQualSbonds(std::vector<std::tuple<std::pair<int,int>,std::pair<int,int>,double> > &contactPairList,double weightCoeff){
  double netLennyJ=0.0;
  if(contactPairList.size()>1){
  for(int l=0;l<contactPairList.size();l++){
    std::tuple<std::pair<int,int>,std::pair<int,int>,double> tple = contactPairList[l];
    double conDist = getPairDistance(std::get<0>(tple),std::get<1>(tple));
    netLennyJ = netLennyJ +lennardJones(std::get<2>(tple),conDist,contactPairList.size(),weightCoeff);
  }
  }
  return netLennyJ;
}

double ktlMolecule::getFitQualSbondsNonSpecific(std::vector<std::pair<std::pair<int,int>,double > > &contactPairList,double weightCoeff){
  double netLennyJ=0.0;
  for(int l=0;l<contactPairList.size();l++){
     double minFit=10000.0;
    for(int m=0;m<contactPairList.size();m++){
      if(m!= l){
         double conDist = getPairDistance(contactPairList[l].first,contactPairList[m].first);
         double lj = lennardJones(contactPairList[l].second,conDist,contactPairList.size(),weightCoeff);
         if(lj<minFit){
           minFit=lj;
         }
        }
    }
    netLennyJ = netLennyJ + minFit;
  }
  return netLennyJ;
}

std::vector<double> ktlMolecule::getExactDistSet(std::vector<std::vector<point> > &coordSet){
  std::vector<double> distSet;
  for(int i=0;i<coordSet.size();i++){
    if(coordSet[i].size()>0){
      for(int j=i;j<coordSet.size();j++){
        if(coordSet[j].size()>0){
        if(j==i){
	          for(int k=0;k<coordSet[i].size()-1;k++){
	            for(int l=k+1;l<coordSet[i].size();l++){
               double d =coordSet[i][k].eDist(coordSet[i][l]);
               distSet.push_back(d);
	            }
    	      }
      }else{
	      for(int k=0;k<coordSet[i].size();k++){
	          for(int l=0;l<coordSet[j].size();l++){
	            double d =coordSet[i][k].eDist(coordSet[j][l]);
	            distSet.push_back(d);
	          }
	        }  
       }
      }
      }
    }
  }
  return distSet;
}


std::vector<double> ktlMolecule::solMolDists(std::vector<std::vector<point> > &pts1){
  std::vector<double> distSet;
  for(int i=0;i<pts1.size();i++){
    for(int j=0;j<pts1[i].size();j++){
      for(int k=0;k<coords.size();k++){
        for(int l=0;l<coords[k].size();l++){
          double d =pts1[i][j].eDist(coords[k][l]);
          distSet.push_back(d);
        }
      }
    }
  }
  return distSet;
}

double  ktlMolecule::allDists(std::vector<point> &ptSet,std::vector<double> &distSet){
  double maxDist=0.0;
  for(int i=0;i<ptSet.size()-1;i++){
    for(int j=i+1;j<ptSet.size();j++){
      double dv =ptSet[i].eDist(ptSet[j]);
      distSet.push_back(dv);
      if(dv>maxDist){
	maxDist=dv;
      }
    }
  }
  return maxDist;
}

std::vector<double> ktlMolecule::getApproxDistSet(double &cutoff){
  std::vector<double> distList;
  // get the mean cloud points
  std::vector<point> centreOfMasses;
  std::vector<double> distSet;
  std::vector<double> maxDists;
  for(int i=0;i<coords.size();i++){
    point meanCoord(0.0,0.0,0.0);
    for(int j=0;j<coords[i].size();j++){
      meanCoord = meanCoord + coords[i][j];
    }
    double maxDist = allDists(coords[i],distSet);
    maxDists.push_back(maxDist);
    meanCoord = meanCoord*(1.0/coords[i].size());
    centreOfMasses.push_back(meanCoord);
  }
  std::vector<double> comDists;
  std::vector<std::pair<int,int> > prs;
  for(int i=0;i<centreOfMasses.size()-1;i++){
    for(int j=i+1;j<centreOfMasses.size();j++){
      double comVal = centreOfMasses[i].eDist(centreOfMasses[j]);
      comDists.push_back(comVal);
      std::pair<int,int> pr;
      pr.first = i;pr.second = j;
      prs.push_back(pr);
    }
  }
  std::cout<<"startvls2\n";
  std::default_random_engine generator;
  for(int i=0;i<comDists.size();i++){
    double lp1 = maxDists[prs[i].first]/3.0;
    double lp2 = maxDists[prs[i].second]/3.0;
    double a1,b1,a2,b2;
    //lp1=7.0;lp2=7.0;
    if(comDists[i]>cutoff){
      a1= -lp1;b1=lp1;a2=-lp2;b2=lp2; 
    }else{
      a1= 0.0;b1=3.0*lp1;a2=0.0;b2=3.0*lp2;
      //a1= -lp1;b1=lp1;a2=-lp2;b2=lp2;
    }
    std::uniform_real_distribution<double> realDistribution1(a1,b1);
    std::uniform_real_distribution<double> realDistribution2(a2,b2);
    for(int j=0;j<2*coords[prs[i].first].size()*coords[prs[i].second].size();j++){
      double extraDist1=realDistribution1(generator);
      double extraDist2=realDistribution2(generator);
      double d = comDists[i] + extraDist1 + extraDist2;
      distSet.push_back(d);
    }
  }
  return distSet;
}


std::vector<double> ktlMolecule::getApproxDistSetFixed(){
  std::vector<double> distList;
  // get the mean cloud points
  std::vector<point> centreOfMasses;
  std::vector<double> distSet;
  std::vector<double> maxDists;
  for(int i=0;i<coords.size();i++){
    point meanCoord(0.0,0.0,0.0);
    for(int j=0;j<coords[i].size();j++){
      meanCoord = meanCoord + coords[i][j];
    }
    double maxDist = allDists(coords[i],distSet);
    maxDists.push_back(maxDist);
    meanCoord = meanCoord*(1.0/coords[i].size());
    centreOfMasses.push_back(meanCoord);
  }
  std::vector<double> comDists;
  std::vector<std::pair<int,int> > prs;
  for(int i=0;i<centreOfMasses.size()-1;i++){
    for(int j=i+1;j<centreOfMasses.size();j++){
      double comVal = centreOfMasses[i].eDist(centreOfMasses[j]);
      comDists.push_back(comVal);
      std::pair<int,int> pr;
      pr.first = i;pr.second = j;
      prs.push_back(pr);
    }
  }
  std::cout<<"startvls2\n";
  std::default_random_engine generator;
  for(int i=0;i<comDists.size();i++){
    //double lp1 = maxDists[prs[i].first]/divFac;
    //double lp2 = maxDists[prs[i].second]/divFac;
    double a1,b1,a2,b2;
    double lp1=7.0;double lp2=7.0;
    if(comDists[i]>5.0){
      a1= -lp1;b1=lp1;a2=-lp2;b2=lp2; 
    }else{
      a1= 0.0;b1=lp1;a2=0.0;b2=lp2;
      //a1= -lp1;b1=lp1;a2=-lp2;b2=lp2;
    }
    std::uniform_real_distribution<double> realDistribution1(a1,b1);
    std::uniform_real_distribution<double> realDistribution2(a2,b2);
    for(int j=0;j<2*coords[prs[i].first].size()*coords[prs[i].second].size();j++){
      double extraDist1=realDistribution1(generator);

      double extraDist2=realDistribution2(generator);
      double d = comDists[i] + extraDist1 + extraDist2;
      distSet.push_back(d);
    }
  }
  return distSet;
}


void ktlMolecule::loadContactPredictions(const char* contactloc){
  std::ifstream cpfile;
  cpfile.open(contactloc);
  if(cpfile.is_open()){
     std::string output;
     while(!cpfile.eof()){
       int ind1,ind2,ind3,ind4;double distance;
       std::getline(cpfile,output);
       std::stringstream ss(output);
       ss>>ind1;
       ss.ignore();
       ss>>ind2;
       std::pair<int,int> pr1;
       pr1.first=ind1;pr1.second=ind2;
       ss.ignore();
       ss>>ind3;
       ss.ignore();
       ss>>ind4;
       ss.ignore();
       std::pair<int,int> pr2;
       pr2.first=ind3;pr2.second=ind4;
       ss>>distance;
       std::tuple<std::pair<int,int>,std::pair<int,int>,double> tp;
       std::get<0>(tp) = pr1;std::get<1>(tp) = pr2;std::get<2>(tp) = distance;
       contactPairList.push_back(tp);
     }
     int cpls =contactPairList.size();
  }	
}

double ktlMolecule::getLennardJonesContact(){
  double ljval =0.0;
  for(int i=0;i<contactPairList.size();i++){
    std::tuple<std::pair<int,int>,std::pair<int,int>,double> tp = contactPairList[i];
    std::pair<int,int> pr1 =  std::get<0>(tp);
    std::pair<int,int> pr2 =  std::get<1>(tp);
    point cd1 = coords[pr1.first][pr1.second];
    point cd2 = coords[pr2.first][pr2.second];
    double prWiseDist = cd1.eDist(cd2);
    double drat = prWiseDist/std::get<2>(tp);
    double indrat  = 1.0/(std::sqrt(2.0)*drat);
    ljval =ljval+ (4.0*(indrat*indrat*indrat*indrat -indrat*indrat)+1.0);
  }
  if(contactPairList.size()>0){
    return ljval;
  }else{
    return ljval;
  }
}
void ktlMolecule::loadFixedSections(const char* fixedsecloc){
  std::ifstream fsfile;
  fsfile.open(fixedsecloc);
  if(fsfile.is_open()){
    std::string output;
     while(!fsfile.eof()){
       std::getline(fsfile,output);
       std::stringstream ss(output);
       int section;
       ss>>section;
       unchangedSections.push_back(section);
     }
  }
}
