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


void Figures::Figure_2()
{

// Set up the data file, can either call COMPOSITE data or CF2 data

Get_data data;
data.get_COMPOSITE();  // alternative:  data.get_CF2();


int nb=11; // define number of shells we want

// Create an hvar object
Hvar hvar1(data,371,264.14,48.26);
// set the direction for this object
hvar1.set_direction(0,10,10);
// run the calculation, this has arguments of (inner cutoff, number of shells, data object)
hvar1.calc_hvar(0,nb);
hvar1.calc_hvar_add(6.25, nb, nb);

ofstream myfile;
myfile.open ("../Figures/data/Figure_2.txt");
double fp;
int status=0;
for (int l=1;l<361;l++)
{

  for (int b=89;b>-91;b--)
  {
  
    
    Hvar hvar2(data,371,264.14,48.26);
    hvar2.set_direction(635,l,b);
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
    
    
    myfile << l << "  " << b << "  " << fp << " " << endl;
    }
    status=(float(l)/360)*100;
    cout<< "\r" << "Generating data . . . " << status << "% complete ";
    std::cout << std::flush;
    //cout<< l << endl;
}

cout << "\n";

myfile.close();

system("python ../Figures/Figure_2.py");

}


