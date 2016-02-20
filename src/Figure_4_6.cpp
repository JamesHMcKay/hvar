#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>

#include <fstream>

#include "boost_offset.hpp"
#include "hvar.hpp"
#include "get_data.hpp"
#include "scanner.hpp"

#include "Figures.hpp"

using namespace std;


void Figures::Figure_4_6()
{

// Set up the data file, can either call COMPOSITE data or CF2 data

Get_data data;
data.get_COMPOSITE();  // alternative:  data.get_CF2();


ofstream myfile;
myfile.open ("../Figures/data/Figure_4.txt");
ofstream myfile2;
myfile2.open ("../Figures/data/Figure_6_b.txt");

ofstream myfile_lnB;
myfile_lnB.open ("../Figures/data/Figure_4_lnB.txt");



Hvar hvarMV(data);
Scanner scan2(hvarMV);
int nb=11;


hvarMV.set_direction(740.6,59.3,16.6);
hvarMV.calc_hvar(6.25,nb);

double pMV;
pMV=hvarMV.prob;



Scanner scan(data);//371,264.14,48.26)




double v_magnitude=0;
double l=0,b=0;
int n=30;
int status;
for (int i=0;i<n;i++)
{
v_magnitude=(float(i)/n)*1600;
scan.find_min_chi_a(v_magnitude,l,b,6.25,11);
l=scan.m_l_conv;
b=scan.m_b_conv;
myfile << v_magnitude << " " << l <<  "  " <<  b <<  "  " << scan.m_min_chi_a << endl;

myfile2 << v_magnitude << " " << l <<  "  " <<  b <<  "  " << scan.m_min_chi_b << endl;

Hvar hvar1(data);
hvar1.set_direction(v_magnitude,l,b);
hvar1.calc_hvar(6.25,nb);

myfile_lnB << v_magnitude << " " << l <<  "  " <<  b <<  "  " << log(pMV/hvar1.prob) << endl;


status=(float(i)/float(n))*50;
cout<< "\r" << "Generating data . . . " << status << "% complete ";
std::cout << std::flush;
}
myfile.close();
myfile2.close();
myfile_lnB.close();


ofstream myfile3;
myfile3.open ("../Figures/data/Figure_5.txt");
ofstream myfile4;
myfile4.open ("../Figures/data/Figure_6_a.txt");

ofstream myfile_lnB_2;
myfile_lnB_2.open ("../Figures/data/Figure_5_lnB.txt");

for (int i=0;i<n;i++)
{
v_magnitude=(float(i)/n)*1200;
scan.find_min_chi_b(v_magnitude,l,b,6.25,11);
l=scan.m_l_conv;
b=scan.m_b_conv;
myfile3 << v_magnitude << " " << l <<  "  " <<  b <<  "  " << scan.m_min_chi_b << endl;
myfile4 << v_magnitude << " " << l <<  "  " <<  b <<  "  " << scan.m_min_chi_a << endl;

Hvar hvar1(data);
hvar1.set_direction(v_magnitude,l,b);
hvar1.calc_hvar(6.25,nb);

myfile_lnB_2 << v_magnitude << " " << l <<  "  " <<  b <<  "  " << log(pMV/hvar1.prob) << endl;



status=(float(i)/float(n))*50+50;
cout<< "\r" << "Generating data . . . " << status << "% complete ";
std::cout << std::flush;
}
myfile.close();
myfile2.close();
myfile_lnB_2.close();


cout << "\n";

system("python ../Figures/Figure_4_6.py");

}
