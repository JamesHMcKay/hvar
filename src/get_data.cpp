#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>

#include <fstream>

#include "boost_offset.hpp"
#include "get_data.hpp"

using namespace std;



void Get_data::get_COMPOSITE()
{



std::ifstream input("../Data/COMPOSITEn-survey-dsrt.dat");
int n=0;
std::string line;
while(getline(input, line)) {
      if (!line.length() || line[0] == '#')
         continue;
      std::istringstream iss(line);
      n=n+1;
   }
  
np=n;


n=0;
std::ifstream input2("../Data/COMPOSITEn-survey-dsrt.dat");
std::string line2;
while(getline(input2, line2)) {
    if (!line2.length() || line2[0] == '#')
       continue;
    std::istringstream iss2(line2);
  
  
  CMBcz.push_back(1);
  d.push_back(1);
  CMBcz.push_back(1);
  CMBvpec.push_back(1);
  sig.push_back(1);
  gl.push_back(1);
  gb.push_back(1);
  iss2>>CMBcz[n] >> d[n] >> CMBvpec[n] >> sig[n] >> gl[n] >> gb[n];
    n=n+1;
 }
 n=n-1;
}


void Get_data::get_CF2()
{
std::ifstream input("../Data/EDD_CF2_data.txt");
int n=0;
std::string line;
while(getline(input, line)) {
      if (!line.length() || line[0] == '#')
         continue;
      std::istringstream iss(line);
      n=n+1;
   }
  
np=n;

double H0=100;
n=0;
double filler;
std::ifstream input2("../Data/EDD_CF2_data.txt");
std::string line2;
while(getline(input2, line2)) {
    if (!line2.length() || line2[0] == '#')
       continue;
    std::istringstream iss2(line2);
  
  
  d.push_back(1);
  CMBcz.push_back(1);
  sig.push_back(1);
  gl.push_back(1);
  gb.push_back(1);
  iss2 >> filler >> d[n] >> sig[n] >> gl[n] >> gb[n] >> CMBcz[n];
  
  d[n]=d[n]*0.745;
  sig[n]=sig[n]*d[n]*H0;
  
    n=n+1;
 }
 n=n-1;
}


