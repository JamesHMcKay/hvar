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



void LG_CMB(Get_data data)
{


// LG compared to CMB reference frame, boost offset 11 shells //

//////////////////////////////////////////////////////////////////////
int nb=11;
//hvar1.calc_hvar_add(6.25,nb,nb);  // add the primed shells to the Hvar object, comment out if not required
// calculate the boost offset between the two hvar objects
//// Create an hvar object
Hvar hvar1(data,371,264.14,48.26);
//// set the direction for this object
hvar1.set_direction(0,0,0);  // can alternatively set direction, wrt to LG if no option given above, or with that specified above (see documentation for explanation)
//// run the calculation, this has arguments of (inner cutoff, number of shells)
hvar1.calc_hvar(6.25,nb);
Hvar hvar2(data);  // no direction specified, default is LG
hvar2.set_direction(0,0,0);//170.869911979 ,  59.1469315251 ,  -4.30994042038 );   //(635,269,28);
hvar2.calc_hvar(6.25,nb);
Boost_offset boost_offset(1);  // enter 1 here to use the model independent calcuation, otherwise set to 0
boost_offset.calc_boost_offset(hvar1,hvar2);
double upper_A=boost_offset.upper_A;
double lower_A=boost_offset.lower_A;
cout<< " ------------------------------------------- " << endl;
cout<< "LG compared to CMB reference frame, 11 shells" << endl;
cout<< " ------------------------------------------- " << endl;
cout<< "Using only primed shells: " << endl;
cout << "calculated b is "<< boost_offset.Get_b() << " +/- " << boost_offset.Get_sig_b() << endl;//  " +/- " << boost_offset.Get_sig_b_systematic() << endl;
cout << "Boost magnitude is "<< boost_offset.Get_v() << " + " << pow(2*(upper_A),0.5)- boost_offset.Get_v() << " - " << boost_offset.Get_v()-pow(2*(lower_A),0.5) << endl;

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//hvar1.calc_hvar_add(6.25,nb,nb);  // add the primed shells to the Hvar object, comment out if not required
// calculate the boost offset between the two hvar objects
//// Create an hvar object
Hvar hvar3(data,371,264.14,48.26);
//// set the direction for this object
hvar3.set_direction(0,0,0);  // can alternatively set direction, wrt to LG if no option given above, or with that specified above (see documentation for explanation)
//// run the calculation, this has arguments of (inner cutoff, number of shells)
hvar3.calc_hvar(0,nb);
Hvar hvar4(data);  // no direction specified, default is LG
hvar4.set_direction(0,0,0);//170.869911979 ,  59.1469315251 ,  -4.30994042038 );   //(635,269,28);
hvar4.calc_hvar(0,nb);

boost_offset.calc_boost_offset(hvar3,hvar4);
upper_A=boost_offset.upper_A;
lower_A=boost_offset.lower_A;
cout<< " ------------------------------------------- " << endl;
cout<< "Using only unprimed shells: " << endl;
cout << "calculated b is "<< boost_offset.Get_b() << " +/- " << boost_offset.Get_sig_b() << endl;//  " +/- " << boost_offset.Get_sig_b_systematic() << endl;
cout << "Boost magnitude is "<< boost_offset.Get_v() << " + " << pow(2*(upper_A),0.5)- boost_offset.Get_v() << " - " << boost_offset.Get_v()-pow(2*(lower_A),0.5) << endl;



//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
hvar3.calc_hvar_add(6.25,nb,nb);
hvar4.calc_hvar_add(6.25,nb,nb);
//// run the calculation, this has arguments of (inner cutoff, number of shells)
Boost_offset boost_offset_2(1);  // enter 1 here to use the model independent calcuation, otherwise set to 0
boost_offset_2.calc_boost_offset(hvar3,hvar4);
upper_A=boost_offset.upper_A;
lower_A=boost_offset.lower_A;
cout<< " ------------------------------------------- " << endl;
cout<< "Using both primed and unprimed shells: " << endl;
cout << "calculated b is "<< boost_offset_2.Get_b() << " +/- " << boost_offset_2.Get_sig_b() << endl;//  " +/- " << boost_offset.Get_sig_b_systematic() << endl;
cout << "Boost magnitude is "<< boost_offset_2.Get_v() << " + " << pow(2*(upper_A),0.5)- boost_offset_2.Get_v() << " - " << boost_offset_2.Get_v()-pow(2*(lower_A),0.5) << endl;


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//hvar1.calc_hvar_add(6.25,nb,nb);  // add the primed shells to the Hvar object, comment out if not required
// calculate the boost offset between the two hvar objects
//// Create an hvar object
Hvar hvar5(data,371,264.14,48.26);
//// set the direction for this object
hvar5.set_direction(0,0,0);  // can alternatively set direction, wrt to LG if no option given above, or with that specified above (see documentation for explanation)
//// run the calculation, this has arguments of (inner cutoff, number of shells)
hvar5.calc_hvar(12.5,nb-1);
Hvar hvar6(data);  // no direction specified, default is LG
hvar6.set_direction(0,0,0);//170.869911979 ,  59.1469315251 ,  -4.30994042038 );   //(635,269,28);
hvar6.calc_hvar(12.5,nb-1);

boost_offset.calc_boost_offset(hvar5,hvar6);
upper_A=boost_offset.upper_A;
lower_A=boost_offset.lower_A;
cout<< " ------------------------------------------- " << endl;
cout<< "Using unprimed shells with inner most shell excluded: " << endl;
cout << "calculated b is "<< boost_offset.Get_b() << " +/- " << boost_offset.Get_sig_b() << endl;//  " +/- " << boost_offset.Get_sig_b_systematic() << endl;
cout << "Boost magnitude is "<< boost_offset.Get_v() << " + " << pow(2*(upper_A),0.5)- boost_offset.Get_v() << " - " << boost_offset.Get_v()-pow(2*(lower_A),0.5) << endl;




//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//hvar1.calc_hvar_add(6.25,nb,nb);  // add the primed shells to the Hvar object, comment out if not required
// calculate the boost offset between the two hvar objects
//// Create an hvar object
Hvar hvar7(data,371,264.14,48.26);
//// set the direction for this object
hvar7.set_direction(0,0,0);  // can alternatively set direction, wrt to LG if no option given above, or with that specified above (see documentation for explanation)
//// run the calculation, this has arguments of (inner cutoff, number of shells)
hvar7.calc_hvar(6.25,nb);
Hvar hvar8(data);  // no direction specified, default is LG
hvar8.set_direction(0,0,0);//170.869911979 ,  59.1469315251 ,  -4.30994042038 );   //(635,269,28);
hvar8.calc_hvar(6.25,nb);

boost_offset.calc_boost_offset_systematic(hvar7,hvar8,nb);
upper_A=boost_offset.upper_A;
lower_A=boost_offset.lower_A;
cout<< " ------------------------------------------- " << endl;
cout<< "With systematic uncertainty (inner cutoff varied from 2-6.25 Mpc/h)" << endl;
cout << "calculated b is "<< boost_offset.Get_b() << " +/- " << boost_offset.Get_sig_b() << " +/- " << boost_offset.Get_sig_b_systematic() << endl;
cout << "Boost magnitude is "<< boost_offset.m_v << endl;//" + " << pow(2*(upper_A),0.5)- boost_offset.Get_v() << " - " << boost_offset.Get_v()-pow(2*(lower_A),0.5) << endl;

cout<< " ------------------------------------------- " << endl;
cout<< " ------------------------------------------- " << endl;


}

void calc_main_results()
{
int nb=11;
Get_data data_COMPOSITE;
data_COMPOSITE.get_COMPOSITE();
Get_data data_CF2;
data_CF2.get_CF2();
cout<< " ------------------------------------------- " << endl;
cout<< " Using COMPOSITE data " << endl;
cout<< " ------------------------------------------- " << endl;
LG_CMB(data_COMPOSITE);


cout<< " Using Cosmicflows-2 data " << endl;
cout<< " ------------------------------------------- " << endl;

LG_CMB(data_CF2);


cout<< " Using COMPOSITE data, best fit frame of reference from LG " << endl;
cout<< " ------------------------------------------- " << endl;


Boost_offset boost_offset(1);
Hvar hvar7(data_COMPOSITE);
//// set the direction for this object
hvar7.set_direction(0,0,0);  // can alternatively set direction, wrt to LG if no option given above, or with that specified above (see documentation for explanation)
//// run the calculation, this has arguments of (inner cutoff, number of shells)
hvar7.calc_hvar(12.5,nb);
Hvar hvar8(data_COMPOSITE);  // no direction specified, default is LG
hvar8.set_direction(170.869911979 ,  59.1469315251 ,  -4.30994042038 );   //(635,269,28);
hvar8.calc_hvar(12.5,nb);
boost_offset.calc_boost_offset_systematic(hvar7,hvar8,nb);
double upper_A=boost_offset.upper_A;
double lower_A=boost_offset.lower_A;
cout<< " ------------------------------------------- " << endl;
cout<< "With systematic uncertainty (inner cutoff varied from 2-6.25 Mpc/h)" << endl;
cout << "calculated b is "<< boost_offset.Get_b() << " +/- " << boost_offset.Get_sig_b() << " +/- " << boost_offset.Get_sig_b_systematic() << endl;
cout << "Boost magnitude is "<< boost_offset.Get_v() << " + " << pow(2*(upper_A),0.5)- boost_offset.Get_v() << " - " << boost_offset.Get_v()-pow(2*(lower_A),0.5) << endl;



//// Create an hvar object
Hvar hvar3(data_COMPOSITE);
//// set the direction for this object
hvar3.set_direction(0,0,0);  // can alternatively set direction, wrt to LG if no option given above, or with that specified above (see documentation for explanation)
//// run the calculation, this has arguments of (inner cutoff, number of shells)
hvar3.calc_hvar(0,nb);
Hvar hvar4(data_COMPOSITE);  // no direction specified, default is LG
hvar4.set_direction(170.869911979 ,  59.1469315251 ,  -4.30994042038);//170.869911979 ,  59.1469315251 ,  -4.30994042038 );   //(635,269,28);
hvar4.calc_hvar(0,nb);
boost_offset.calc_boost_offset(hvar3,hvar4);
upper_A=boost_offset.upper_A;
lower_A=boost_offset.lower_A;
cout<< " ------------------------------------------- " << endl;
cout<< "Using only unprimed shells: " << endl;
cout << "calculated b is "<< boost_offset.Get_b() << " +/- " << boost_offset.Get_sig_b() << endl;//  " +/- " << boost_offset.Get_sig_b_systematic() << endl;
cout << "Boost magnitude is "<< boost_offset.Get_v() << " + " << pow(2*(upper_A),0.5)- boost_offset.Get_v() << " - " << boost_offset.Get_v()-pow(2*(lower_A),0.5) << endl;



//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
hvar3.calc_hvar_add(6.25,nb,nb);
hvar4.calc_hvar_add(6.25,nb,nb);
//// run the calculation, this has arguments of (inner cutoff, number of shells)
Boost_offset boost_offset_2(1);  // enter 1 here to use the model independent calcuation, otherwise set to 0
boost_offset_2.calc_boost_offset(hvar3,hvar4);
upper_A=boost_offset.upper_A;
lower_A=boost_offset.lower_A;
cout<< " ------------------------------------------- " << endl;
cout<< "Using both primed and unprimed shells: " << endl;
cout << "calculated b is "<< boost_offset_2.Get_b() << " +/- " << boost_offset_2.Get_sig_b() << endl;//  " +/- " << boost_offset.Get_sig_b_systematic() << endl;
cout << "Boost magnitude is "<< boost_offset_2.m_v << endl;// << " + " << pow(2*(upper_A),0.5)- boost_offset_2.Get_v() << " - " << boost_offset_2.Get_v()-pow(2*(lower_A),0.5) << endl;


cout<< " ------------------------------------------- " << endl;
cout<< " ------------------------------------------- " << endl;
cout<< " ln(B) values for LG versus CMB " << endl;
}



int main()
{
int nb=11;
Get_data data_COMPOSITE;
data_COMPOSITE.get_COMPOSITE();
Get_data data_CF2;
data_CF2.get_CF2();

calc_main_results();


Hvar hvarCMB(data_COMPOSITE,371,264.14,48.26);
//// set the direction for this object
hvarCMB.set_direction(0,0,0);  // can alternatively set direction, wrt to LG if no option given above, or with that specified above (see documentation for explanation)
//// run the calculation, this has arguments of (inner cutoff, number of shells)
hvarCMB.calc_hvar(6.25,nb);
Hvar hvarLG(data_COMPOSITE);  // no direction specified, default is LG
hvarLG.set_direction(0,0,0);//170.869911979 ,  59.1469315251 ,  -4.30994042038 );   //(635,269,28);
hvarLG.calc_hvar(6.25,nb);
//

Hvar hvarMV(data_COMPOSITE);  // no direction specified, default is LG
hvarMV.set_direction(740,59.3,16.6);//170.869911979 ,  59.1469315251 ,  -4.30994042038 );   //(635,269,28);
hvarMV.calc_hvar(6.25,nb);

Hvar hvarB(data_COMPOSITE);  // no direction specified, default is LG
hvarB.set_direction(170.869911979 ,  59.1469315251 ,  -4.30994042038 ); // B for the best systematic boost frame
hvarB.calc_hvar(6.25,nb);



double pCMB;
double pLG;
double pMV;
double pB;
pLG=hvarLG.prob;
pCMB=hvarCMB.prob;
pMV=hvarMV.prob;
pB=hvarB.prob;

cout<< "ln B = (pLG/pCMB) " << log(pLG/pCMB) << endl;

cout<< "ln B = (pMV/pCMB) " << log(pMV/pCMB) << endl;

cout<< "ln B = (pMV/pLG) " << log(pMV/pLG) << endl;

cout<< "ln B = (pMV/pB) " << log(pMV/pB) << endl;
cout<< "chi2_b for MV " << hvarMV.m_chi2_a <<endl;

cout<< "chi2_b for LG " << hvarLG.m_chi2_a <<endl;

cout<< "chi2_b for CMB " << hvarCMB.m_chi2_a <<endl;




}



// how to print out chi^2 if required
//cout<< "chi2_a for LG " << hvar1.m_chi2_a <<endl;
//cout<< "chi2_a for CMB " << hvar2.m_chi2_a <<endl;  // print out chi^2 if required

// example of how to call the Hubble parameter and average shell radius from a hvar object
//std::vector<double> H;
//std::vector<double> r;
//H=hvar1.m_H;
//r=hvar1.m_av_r2;

