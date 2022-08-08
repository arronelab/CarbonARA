#include "localWrithe.h"
#include <time.h>
#include <string.h>


int main( int argc, const char* argv[] )
{
	bool check=false;
	std::ifstream myfile;
	std::ifstream myfile2;
 	myfile.open(argv[1]);
	myfile2.open(argv[2]);
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
	std::string output2;
	std::vector<point> points2;
	if (myfile2.is_open()) { 
	  while(std::getline(myfile2,output2)){
	    points2.push_back(point(output2));
	  }
	}else{
	  std::cout<<"Curve data file 2 failed to open";
	}
	if(points2.size()>3){
	  check =true;
	}
	myfile2.close();
	if(check){
	  localWrithe lw;
	  std::cout<<lw.DIClosedLk(points,points2)<<"\n";
	}else{
	  std::cout<<"must supply a file containing the curve data";
	}
	return 0;
}
       
