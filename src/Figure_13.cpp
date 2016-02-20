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


void Figures::Figure_13()
{

// Set up the data file, can either call COMPOSITE data or CF2 data

Get_data data;
data.get_CF2();  // alternative:  data.get_CF2();


int nb=11; // define number of shells we want

// Create an hvar object
Hvar hvar1(data);
// set the direction for this object
hvar1.set_direction(0,10,10);
// run the calculation, this has arguments of (inner cutoff, number of shells, data object)
hvar1.calc_hvar(6.25,nb);
hvar1.calc_hvar_add(0,nb,nb);

// repeat for a second direction
Hvar hvar2(data,371,264.14,48.26);
hvar2.set_direction(0,0,0);
hvar2.calc_hvar(6.25,nb);
hvar2.calc_hvar_add(0,nb,nb);



// example of how to call the Hubble parameter and average shell radius from a hvar object
std::vector<double> H1;
std::vector<double> r2;
std::vector<double> H2;
std::vector<double> sig_H2;
std::vector<double> sig_H1;
std::vector<double> sig_r2;
std::vector<double> sig_Del_H(nb*2);
sig_H1=hvar1.m_sig_H;
sig_H2=hvar2.m_sig_H;
H1=hvar1.m_H;
H2=hvar2.m_H;
r2=hvar1.m_av_r2;
sig_r2=hvar1.m_sig_av_r2;

for (int i=0;i<nb*2;i++)
{
sig_Del_H[i]=pow( pow(sig_H1[i],2)+pow(sig_H2[i],2) , 0.5);
}
ofstream myfile;
myfile.open ("../Figures/data/Figure_13.txt");
for (int i=0;i<nb*2;i++)
{
myfile << r2[i] << "  " << H1[i]-H2[i] << "  " << sig_r2[i] << " " << sig_Del_H[i] <<  endl;

}
myfile.close();
// calculate the boost offset between the two hvar objects
ofstream myfile2;
myfile2.open ("../Figures/data/Figure_13_details.txt");
Boost_offset boost_offset;
boost_offset.calc_boost_offset(hvar1,hvar2);
myfile2 <<"  " << boost_offset.Get_b() << " " << boost_offset.Get_A() <<  endl;

// now repeat for unprimed shells
Hvar hvar3(data);
hvar3.set_direction(0,10,10);
hvar3.calc_hvar(6.25,nb);
Hvar hvar4(data,371,264.14,48.26);
hvar4.set_direction(0,0,0);
hvar4.calc_hvar(6.25,nb);
//Boost_offset boost_offset;
boost_offset.calc_boost_offset(hvar3,hvar4);
myfile2 <<"  " << boost_offset.Get_b() << " " << boost_offset.Get_A() <<  endl;
// now repeat for primed shells
Hvar hvar5(data);
hvar5.set_direction(0,10,10);
hvar5.calc_hvar(0,nb);
Hvar hvar6(data,371,264.14,48.26);
hvar6.set_direction(0,0,0);
hvar6.calc_hvar(0,nb);
//Boost_offset boost_offset;
boost_offset.calc_boost_offset(hvar5,hvar6);
myfile2 <<"  " << boost_offset.Get_b() << " " << boost_offset.Get_A() <<  endl;




myfile2.close();

system("python ../Figures/Figure_13.py");

}


