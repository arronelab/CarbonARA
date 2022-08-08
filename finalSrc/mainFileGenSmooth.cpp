#include "localWrithe.h"
#include <time.h>


int main( int argc, const char* argv[] )
{
	if(argc ==2){
	bool check=false;
	std::ifstream myfile;
 	myfile.open(argv[1]);
	std::string output;
	std::vector<point> points;
	if (myfile.is_open()) { 
	  while(std::getline(myfile,output)){
		points.push_back(point(output));
  		}
	}else{
	std::cout<<"Curve data file failed to open";
	}
	if(points.size()>3){
	    check =true;
	}
	myfile.close();
	if(check){
	localWrithe lw;
	std::vector<double> writheOutput;
	writheOutput = lw.getWritheGenSmooth(points);
	std::cout<<writheOutput[0]<<" "<<writheOutput[1]<<" "<<writheOutput[2]<<" {";
	if(writheOutput.size()>3){
	  for(int i=3;i<writheOutput.size()-1;i++){
	    std::cout<<writheOutput[i]<<" ";
	  }
	  std::cout<<writheOutput[writheOutput.size()-1]<<"}\n";
	}else{
	  std::cout<<"}";
	}
	}
	}else{
	  std::cout<<"must supply a file containing the curve data";
	}
	return 0;
}
