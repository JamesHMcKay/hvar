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


void Figures::Figure_8()
{
int use_H_bar=1;
// Set up the data file, can either call COMPOSITE data or CF2 data

Get_data data;
data.get_COMPOSITE();  // alternative:  data.get_CF2();
//
//
int nb=11; // define number of shells we want
//
//// Create an hvar object
Hvar hvar1(data);//,371,264.14,48.26);
//// set the direction for this object
hvar1.set_direction(0,10,10);
//// run the calculation, this has arguments of (inner cutoff, number of shells, data object)
hvar1.calc_hvar(0,nb);
hvar1.calc_hvar_add(6.25,nb,nb);
//



Scanner scan(hvar1,use_H_bar);

ofstream myfile;
myfile.open ("../Figures/data/Figure_8_all_shells.txt");
double ratio=0;
int status;

const int n=20;
double v_magnitude[n+1];
double l[n+1],b[n+1];
double S;
l[0]=0;
b[0]=0;
int burnin_length=2000;
double sig_l[n+1],sig_b[n+1];
for (int i=1;i<n+1;i++)
{
v_magnitude[i]=(float(i-1)/n)*400+80;
scan.find_best_fit(v_magnitude[i],l[i-1],b[i-1],burnin_length,3000);
Boost_offset boost_offset=scan.Get_boost_object();
l[i]=scan.m_l_conv;
b[i]=scan.m_b_conv;
sig_l[i]=scan.m_l_sig;
sig_b[i]=scan.m_b_sig;
ratio=boost_offset.Get_v()/v_magnitude[i];
S=boost_offset.Get_chi_S();
myfile << v_magnitude[i] << " " << l[i] <<  "  " <<  b[i] <<  "  " << ratio << " " << S << " " << sig_l[i] << " " << sig_b[i] << endl;
status=(float(i-1)/float(n))*100;
cout<< "\r" << "Generating data . . . " << status << "% complete ";
std::cout << std::flush;
system("gnuplot ../Figures/histogram_l.gnu");
system("gnuplot ../Figures/histogram_b.gnu");
burnin_length=1000;
}
myfile.close();



// calculate quantities for primed and unprimed shells only, using same (l,b) values

ofstream myfile_primed;
myfile_primed.open ("../Figures/data/Figure_8_all_primed.txt");
ofstream myfile_unprimed;
myfile_unprimed.open ("../Figures/data/Figure_8_all_unprimed.txt");

Hvar hvar1_primed(data);
hvar1_primed.set_direction(0,10,10);
hvar1_primed.calc_hvar(6.25,nb);

Hvar hvar1_unprimed(data);
hvar1_unprimed.set_direction(0,10,10);
hvar1_unprimed.calc_hvar(0,nb);


for (int i=1;i<n+1;i++)
{

Hvar hvar2_primed(data);//371,264.14,48.26);
hvar2_primed.set_direction(v_magnitude[i],l[i],b[i]);   //(635,269,28);
hvar2_primed.calc_hvar(6.25,nb);

Hvar hvar2_unprimed(data);//371,264.14,48.26);
hvar2_unprimed.set_direction(v_magnitude[i],l[i],b[i]);   //(635,269,28);
hvar2_unprimed.calc_hvar(0,nb);

Boost_offset boost_offset_primed(use_H_bar);

boost_offset_primed.calc_boost_offset(hvar1_primed,hvar2_primed);

Boost_offset boost_offset_unprimed(use_H_bar);

boost_offset_unprimed.calc_boost_offset(hvar1_unprimed,hvar2_unprimed);

ratio=boost_offset_primed.Get_v()/v_magnitude[i];
S=boost_offset_primed.Get_chi_S();
myfile_primed << v_magnitude[i] << " " << l[i] <<  "  " <<  b[i] <<  "  " << ratio << " " << S << " " << sig_l[i] << " " << sig_b[i] << endl;



ratio=boost_offset_unprimed.Get_v()/v_magnitude[i];
S=boost_offset_unprimed.Get_chi_S();
myfile_unprimed << v_magnitude[i] << " " << l[i] <<  "  " <<  b[i] <<  "  " << ratio << " " << S << " " << sig_l[i] << " " << sig_b[i] << endl;



}

myfile_primed.close();

myfile_unprimed.close();



system("python ../Figures/Figure_8.py");


cout << "\n";



}


