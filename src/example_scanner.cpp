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

using namespace std;

int main()
{

// Set up the data file, can either call COMPOSITE data or CF2 data

Get_data data;
data.get_CF2();  // alternative:  data.get_CF2();
//
//
int nb=11; // define number of shells we want
int burnin_length=1000; // define length of burn in
int chain_length=3000; // define number of steps after burn in complete
//
//// Create an hvar object
Hvar hvar1(data); //371,264.14,48.26);
//// set the direction for this object
hvar1.set_direction(0,0,0);
//// fill the Hvar object
hvar1.calc_hvar(0,nb);
hvar1.calc_hvar_add(6.25,nb,nb);



// initialise the Scanner object
Scanner scan(hvar1,0);

scan.find_best_fit(400,0,0,burnin_length,chain_length);  // call the scanner, inputs are (v,l,b)
//Boost_offset boost_offset=scan.Get_boost_object();

// plot the distrbution in the l and b parameter spaces
system("gnuplot ../Figures/histogram_l.gnu");
system("gnuplot ../Figures/histogram_b.gnu");

}


