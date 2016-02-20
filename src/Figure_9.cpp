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


void Figures::Figure_9()
{



int n=0;
std::ifstream input2("../Figures/data/Figure_8_all_shells.txt");
std::string line2;
while(getline(input2, line2)) {
      if (!line2.length() || line2[0] == '#')
         continue;
      std::istringstream iss(line2);
      n=n+1;
   }
 
if (n==0)
{
cout << "Must run Figure_8 first, please compile and run Figure_8 then try again"<<endl;
return;
}

double v_bestfit,l_bestfit,b_bestfit;

system("python ../Figures/Figure_9_a.py");

std::ifstream input("../Figures/data/Figure_8_best_fit.txt");
std::string line;
while(getline(input, line)) {
      if (!line.length() || line[0] == '#')
         continue;
      std::istringstream iss(line);
      iss>> v_bestfit >> l_bestfit >> b_bestfit;
   }


// Set up the data file, can either call COMPOSITE data or CF2 data

Get_data data;
data.get_COMPOSITE();  // alternative:  data.get_CF2();



int nb=11; // define number of shells we want

// Create an hvar object
Hvar hvar1(data);
// set the direction for this object
hvar1.set_direction(0,0,0);
// run the calculation, this has arguments of (inner cutoff, number of shells, data object)
hvar1.calc_hvar(0,nb);
hvar1.calc_hvar_add(6.25,nb,nb);

// repeat for a second direction
Hvar hvar2(data);
hvar2.set_direction(v_bestfit,l_bestfit,b_bestfit);
hvar2.calc_hvar(0,nb);
hvar2.calc_hvar_add(6.25,nb,nb);



// example of how to call the Hubble parameter and average shell radius from a hvar object
std::vector<double> H1;
std::vector<double> r2;
std::vector<double> H2;
std::vector<double> sig_H2;
std::vector<double> sig_H1;
std::vector<double> sig_r2;
std::vector<double> r;

std::vector<double> dH1;
std::vector<double> dH2;
std::vector<double> sig_dH2;
std::vector<double> sig_dH1;

std::vector<double> sig_Del_H(nb*2);
sig_H1=hvar1.m_sig_H;
sig_H2=hvar2.m_sig_H;
H1=hvar1.m_H;
H2=hvar2.m_H;
r2=hvar1.m_av_r2;
r=hvar1.m_av_r;

dH1=hvar1.m_delta_H;
dH2=hvar2.m_delta_H;
sig_dH1=hvar1.m_sig_delta_H;
sig_dH2=hvar2.m_sig_delta_H;

sig_r2=hvar1.m_sig_av_r2;

for (int i=0;i<nb*2;i++)
{
sig_Del_H[i]=pow( pow(sig_H1[i],2)+pow(sig_H2[i],2) , 0.5);
}
ofstream myfile;
myfile.open ("../Figures/data/Figure_9_a.txt");
ofstream myfile2;
myfile2.open ("../Figures/data/Figure_9_b.txt");
ofstream myfile3;
myfile3.open ("../Figures/data/Figure_9_c.txt");



for (int i=0;i<nb*2;i++)
{
myfile << r2[i] << "  " << H1[i]-H2[i] << "  " << sig_r2[i] << " " << sig_Del_H[i] <<  endl;
myfile2 << r[i] << "  " << dH1[i] << " " << sig_dH1[i] <<  endl;
myfile3 << r[i] << "  " << dH2[i] << " " << sig_dH2[i] <<  endl;
}
myfile.close();
myfile2.close();
myfile3.close();
// calculate the boost offset between the two hvar objects
ofstream myfile4;
myfile4.open ("../Figures/data/Figure_9_details.txt");
Boost_offset boost_offset(0); // using 0 since we are ploting, want A model independent not v
boost_offset.calc_boost_offset(hvar1,hvar2);
myfile4 <<"  " << boost_offset.Get_b() << " " << boost_offset.Get_A() <<  endl;

// now repeat for unprimed shells
Hvar hvar3(data);
hvar3.set_direction(0,0,0);
hvar3.calc_hvar(6.25,nb);
Hvar hvar4(data);
hvar4.set_direction(v_bestfit,l_bestfit,b_bestfit);
hvar4.calc_hvar(6.25,nb);
//Boost_offset boost_offset;
boost_offset.calc_boost_offset(hvar3,hvar4);
myfile4 <<"  " << boost_offset.Get_b() << " " << boost_offset.Get_A() <<  endl;
// now repeat for primed shells
Hvar hvar5(data);
hvar5.set_direction(0,0,0);
hvar5.calc_hvar(0,nb);
Hvar hvar6(data);
hvar6.set_direction(v_bestfit,l_bestfit,b_bestfit);
hvar6.calc_hvar(0,nb);
//Boost_offset boost_offset;
boost_offset.calc_boost_offset(hvar5,hvar6);
myfile4 <<"  " << boost_offset.Get_b() << " " << boost_offset.Get_A() <<  endl;




myfile2.close();

system("python ../Figures/Figure_9_b.py");

}


