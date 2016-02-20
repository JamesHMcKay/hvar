#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>

#include <fstream>

#include "boost_offset.hpp"
#include "hvar.hpp"
#include "get_data.hpp"
#include "Figures.hpp"

using namespace std;


void Figures::Figure_3()
{

// Set up the data file, can either call COMPOSITE data or CF2 data

Get_data data;
data.get_COMPOSITE();  // alternative:  data.get_CF2();


int nb=11; // define number of shells we want



ofstream myfile;
myfile.open ("../Figures/data/Figure_3_a.txt");
int status=0;
int df=8;
for (int l=1;l<361;l++)
{

  //for (int b=-89;b<91;b++)
  for (int b=89;b>-91;b--)
  {
  // Create an hvar object
    Hvar hvar1(data);
    hvar1.set_direction(635,l,b);
    hvar1.calc_hvar(6.25,nb);

    
    myfile << l << "  " << b << "  " << hvar1.m_chi2_a/(df) << " " << endl;
    }
    status=(float(l)/360)*50;
    cout<< "\r" << "Generating data . . . " << status << "% complete ";
    std::cout << std::flush;

}

myfile.close();


system("python ../Figures/Figure_3_a.py");


ofstream myfile2;
myfile2.open ("../Figures/data/Figure_3_b.txt");

for (int l=1;l<361;l++)
{

  //for (int b=-89;b<91;b++)
  for (int b=89;b>-91;b--)
  {
  
   // Create an hvar object
    Hvar hvar1(data);
    hvar1.set_direction(740,l,b);
    hvar1.calc_hvar(6.25,nb);

    
    
    
    myfile2 << l << "  " << b << "  " << hvar1.m_chi2_a/(df) << " " << endl;
    }
    status=(float(l)/360)*50+50;
    cout<< "\r" << "Generating data . . . " << status << "% complete ";
    std::cout << std::flush;

}

myfile2.close();

system("python ../Figures/Figure_3_b.py");

}


