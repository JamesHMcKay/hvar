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
#include "scanner.hpp"

using namespace std;


void Figures::Figure_23()
{

// Set up the data file, can either call COMPOSITE data or CF2 data

Get_data data;
data.get_COMPOSITE();  // alternative:  data.get_CF2();


int nb=11; // define number of shells we want

int df=1;

ofstream myfile;
myfile.open ("../Figures/data/Figure_23_a.txt");


ofstream myfile_chi2;
myfile_chi2.open ("../Figures/data/Figure_23_a_chi2.txt");

int status=0;

Hvar hvarMV(data);
Scanner scan(hvarMV);


scan.find_min_chi_a(635,0,0,6.25,11);
double l_MV=scan.m_l_conv;
double b_MV=scan.m_b_conv;


hvarMV.set_direction(635,l_MV,b_MV);
hvarMV.calc_hvar(6.25,nb);

double pMV;
pMV=hvarMV.prob;
//double p;

for (int l=1;l<361;l++)
{

  //for (int b=-89;b<91;b++)
  for (int b=89;b>-91;b--)
  {
  // Create an hvar object
    Hvar hvar1(data);
    hvar1.set_direction(635,l,b);
    hvar1.calc_hvar(6.25,nb);
    
    
    myfile << l << "  " << b << "  " << log(pMV/hvar1.prob) << " " << endl;
    
    myfile_chi2 << l << "  " << b << "  " << hvar1.m_chi2_a/(df) << " " << endl;
    }
    status=(float(l)/360)*50;
    cout<< "\r" << "Generating data . . . " << status << "% complete ";
    std::cout << std::flush;

}

myfile.close();
myfile_chi2.close();


system("python ../Figures/Figure_23_a.py");

scan.find_min_chi_a(740,0,0,6.25,11);
l_MV=scan.m_l_conv;
b_MV=scan.m_b_conv;


hvarMV.set_direction(635,l_MV,b_MV);
hvarMV.calc_hvar(6.25,nb);

pMV=hvarMV.prob;

ofstream myfile2;
myfile2.open ("../Figures/data/Figure_23_b.txt");

ofstream myfile2_chi2;
myfile2_chi2.open ("../Figures/data/Figure_23_b_chi2.txt");

for (int l=1;l<361;l++)
{

  //for (int b=-89;b<91;b++)
  for (int b=89;b>-91;b--)
  {
  
   // Create an hvar object
    Hvar hvar1(data);
    hvar1.set_direction(740,l,b);
    hvar1.calc_hvar(6.25,nb);

    
    
    
    myfile2 << l << "  " << b << "  " << log(pMV/hvar1.prob) << " " << endl;
    
    myfile2_chi2 << l << "  " << b << "  " << hvar1.m_chi2_a/(df) << " " << endl;
    }
    status=(float(l)/360)*50+50;
    cout<< "\r" << "Generating data . . . " << status << "% complete ";
    std::cout << std::flush;

}

myfile2.close();
myfile2_chi2.close();

system("python ../Figures/Figure_23_b.py");

}


