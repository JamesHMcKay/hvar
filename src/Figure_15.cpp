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


void Figures::Figure_15()
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
hvar1.calc_hvar(0,nb);
hvar1.calc_hvar_add(6.25, nb, nb);

ofstream myfile;
myfile.open ("../Figures/data/Figure_15_a.txt");
double fp;
int status=0;

double v_1=450,v_2=200;

for (int l=1;l<361;l++)
{
  for (int b=89;b>-91;b--)
  {
    Hvar hvar2(data);
    hvar2.set_direction(v_1,l,b);
    hvar2.calc_hvar(0,nb);
    hvar2.calc_hvar_add(6.25, nb, nb);
    // calculate the boost offset between the two hvar objects
    Boost_offset boost_offset;
    boost_offset.calc_boost_offset(hvar2,hvar1);
    if (boost_offset.Get_sign_A()==0)
    {
    fp=abs(boost_offset.Get_b()+1);
    }
    else
    {
    fp=2-abs(boost_offset.Get_b()+1);
    }
    
    if (fp>2)
    {
    fp=2;
    }
    myfile << l << "  " << b << "  " << fp << " " << endl;
    Hvar hvar3(data);
    hvar3.set_direction(v_1,l,b);
    hvar3.calc_hvar(6.25,nb);
    }
    status=(float(l)/360)*50;
    cout<< "\r" << "Generating data . . . " << status << "% complete ";
    std::cout << std::flush;
}
cout << "\n";
myfile.close();
system("python ../Figures/Figure_15_a.py");


ofstream myfile3;
myfile3.open ("../Figures/data/Figure_15_b.txt");

for (int l=1;l<361;l++)
{
  for (int b=89;b>-91;b--)
  {
    Hvar hvar2(data);
    hvar2.set_direction(v_2,l,b);
    hvar2.calc_hvar(0,nb);
    hvar2.calc_hvar_add(6.25, nb, nb);
    // calculate the boost offset between the two hvar objects
    Boost_offset boost_offset;
    boost_offset.calc_boost_offset(hvar2,hvar1);
    if (boost_offset.Get_sign_A()==0)
    {
    fp=abs(boost_offset.Get_b()+1);
    }
    else
    {
    fp=2-abs(boost_offset.Get_b()+1);
    }
    
    if (fp>2)
    {
    fp=2;
    }
    myfile3 << l << "  " << b << "  " << fp << " " << endl;
    Hvar hvar3(data);
    hvar3.set_direction(v_2,l,b);
    hvar3.calc_hvar(6.25,nb);
    }
    status=(float(l)/360)*50+50;
    cout<< "\r" << "Generating data . . . " << status << "% complete ";
    std::cout << std::flush;
}
cout << "\n";
myfile3.close();
system("python ../Figures/Figure_15_b.py");




}


