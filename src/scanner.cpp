#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include "boost_offset.hpp"
#include "scanner.hpp"

#include <time.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;



void Scanner::find_best_fit(double v,double l,double b,int burnin_length, int chain_length)
{
double minr_a=m_hvar1.m_minr_a;
double minr_b=m_hvar1.m_minr_b;
double l_trial=0,b_trial=0;
double l_sum=0,b_sum=0, tol=10;
double accept_ratio=0, step_size=3;
int burn_in=0;
int nb=m_hvar1.m_nb,two_configurations=0;
int nb_prev=m_hvar1.m_nb_prev,k=0,n=0,accept=0,accepted=0,q=0,bad=0,good=0;
int seed=5,converged=0,stuck=0;
double x,y,z,fp,fp_new=0;
std::vector<double> l_accepted(chain_length),b_accepted(chain_length);
ofstream myfile;
myfile.open ("../Figures/data/Figure_8_test.txt");
ofstream myfile2;
myfile2.open ("../Figures/data/Figure_8_hist_l.txt");
ofstream myfile3;
myfile3.open ("../Figures/data/Figure_8_hist_b.txt");
Ran myran(seed);
int nb_b=0,nb_a=nb;
// check how many shell configurations we need to use
if (nb_prev!=0)
{
nb_b=nb-nb_prev;
nb_a=nb_prev;
two_configurations=1;
}

// function to minimise
Hvar hvar2(m_data,m_v_ref,m_l_ref,m_b_ref);
hvar2.set_direction(v,l,b);
hvar2.calc_hvar(minr_a,nb_a);
if (two_configurations==1)
{
hvar2.calc_hvar_add(minr_b,nb_b,nb_a);
}
Boost_offset boost_offset(m_use_H_bar);
boost_offset.calc_boost_offset(hvar2,m_hvar1);
fp=boost_offset.Get_fp();
double sig_b=pow(boost_offset.Get_sig_b(),0.25);//exp(pow(1-(sqrt(2*70*boost_offset.Get_A())/560),2));

fp=fp*sig_b;
cout << "initial fp = " << fp << endl;


while (converged==0)
{
  
    n=n+1;

    // burn in control
    if (accepted==burnin_length && burn_in==0)
    {
    cout << "\n burn in complete, fp = " << fp << " accepted points = " << q << endl;
    burn_in=1;
    k=0;
    n=0;
    accepted=0;
    }

    // calculate new trial points
    int flag=0; // use this to control the following loop so that we produce a random step for both l and b
    while (flag<2)
    {
        x=myran.doub();
        y=-log(x);
        z=myran.doub();
        if (z<= exp(  -0.5*pow(y-1,2)    )  )
        {
            z=myran.doub();
            if (z<0.5)
            {
            y=-y;
            }
          
            if (fp_new==2*sig_b)
            {
              y=y*500;  // if we are in the wrong side of the sky (fp==2, A<0) make bigger steps
              if (stuck>100)
              {
              y=y/10;
              }
              stuck=stuck+1;
              
              cout<< "\r" << "stuck = " << stuck;
              std::cout << std::flush;
        
            }
            else
            {
            y=y*pow(accept_ratio/0.3,3); // adjust step size to keep acceptance rate appropriate (not get stuck)
            }

            if (flag==0)
            {
              l_trial=fmod(l+y*step_size,360);
              if (l_trial<0)
              {
              l_trial=0;
              }
            }
            else
            {
              b_trial=(b+y*step_size);
              if (b_trial>=90)
              {
               b_trial=90;
              }
              if (b_trial<=-90)
              {
               b_trial=-90;
              }
            }
            flag=flag+1;
        }
    }
  
  
  
    // evaluate fp for this new point
    Hvar hvar2(m_data,m_v_ref,m_l_ref,m_b_ref);
    hvar2.set_direction(v,l_trial,b_trial);
    hvar2.calc_hvar(minr_a,nb_a);
    if (two_configurations==1)
    {
    hvar2.calc_hvar_add(minr_b,nb_b,nb_a);
    }
    Boost_offset boost_offset(m_use_H_bar);
    boost_offset.calc_boost_offset(hvar2,m_hvar1);
    fp_new=boost_offset.Get_fp();
    sig_b=pow(boost_offset.Get_sig_b(),0.25);//exp(pow(1-(sqrt(2*70*boost_offset.Get_A())/560),2));

    fp_new=fp_new*sig_b;
  
  
  
    // accept or reject new point
    accept=0;
    if (fp_new<=fp)
    {
      accept=1;
      good=good+1;
      fp=fp_new;
      l=l_trial;
      b=b_trial;
    }
    else if (abs(y)*myran.doub()<fp/fp_new/tol && fp_new<2*sig_b)
    {
      accept=1;
      bad=bad+1;
      fp=fp_new;
      l=l_trial;
      b=b_trial;
    }
  
    // control the accept ratio
    if (accepted>100)
    {
    accept_ratio=float(accepted)/float(n);
    }



    // store accepted points
    if (accept==1 && fp<2*sig_b)
    {
      
    
      if (burn_in==1)
      {
            myfile2 << l << endl;
            myfile3 << b << endl;
            myfile << q << " " << l << " " << b << " " << fp << endl;
            b_sum=b+b_sum;
            l_sum=l+l_sum;

            l_accepted[k]=l;
            b_accepted[k]=b;
            k=k+1;
        
        
      }

      accepted=accepted+1;
      q=q+1;
      
      cout<< "\r" << "Acceptence ratio = " << accept_ratio << " bad steps = " << bad << " good steps = " << good << " propsed steps = " << n << " fp = " << fp << " step size = " << abs(y) ;
      std::cout << std::flush;

     }
  
  

  
    // converge after set number of accepted points
    if (burn_in==1 && accepted==chain_length)
    {
    converged=1;
    
    //evaluate mean and standard deviation for (l,b)
    m_l_conv=l_sum/k;
    m_b_conv=b_sum/k;
    double var_l=0,var_b=0;
    for (int j=1;j<chain_length;j++)
    {
     var_l=pow(l_accepted[j]-m_l_conv,2)+var_l;
     var_b=pow(b_accepted[j]-m_b_conv,2)+var_b;
    }
    var_b=pow(var_b/k,0.5);
    var_l=pow(var_l/k,0.5);
    
    m_l_sig=var_l;
    m_b_sig=var_b;
      
    // evaluate the boost offset for this averaged point
    Hvar hvar2(m_data,m_v_ref,m_l_ref,m_b_ref);
    hvar2.set_direction(v,m_l_conv,m_b_conv);
    hvar2.calc_hvar(minr_a,nb_a);
    if (two_configurations==1)
    {
    hvar2.calc_hvar_add(minr_b,nb_b,nb_a);
    }
    Boost_offset boost_offset(m_use_H_bar);
    boost_offset.calc_boost_offset(hvar2,m_hvar1);
    minimum_boost=boost_offset;  // store boost offset object for later retrieval
    
    // output results
    cout<< "\n Scanner converged at l = " << m_l_conv << " +/-  " << var_l << " b = "<< m_b_conv << " +/- " << var_b << endl;
    cout << "b = " << boost_offset.Get_b() << " +/-" << boost_offset.Get_sig_b() << endl;
    }
}

myfile.close();
myfile2.close();
myfile3.close();

}



void Scanner::find_min_chi(double v,double l,double b,double minr,int nb,int use_chi2_a)
{


Hvar hvar(m_data,m_v_ref,m_l_ref,m_b_ref);

int seed,converged=0;
double x,y,z;
seed=5;
Ran myran(seed);
double f; // function to minimise
hvar.set_direction(v,l,b);
hvar.calc_hvar(minr,nb);


if (use_chi2_a==1)
{
f=hvar.m_chi2_a;
}
else
{
f=hvar.m_chi2_b;
}


double l_trial=0,b_trial=0,f_new;

//cout << "f = " << f << endl;
int k=0,n=0,burn_in=0;
int accept=0;
double l_sum=0,b_sum=0;
double tol=0.3;
int accepted=0;
double step_size=2;
while (converged==0)
{
    // calculate new trial points
    
        n=n+1;
        if (n>200 && burn_in==0)
        {
        //cout << "\n burn in complete" << endl;
        burn_in=1;
        k=0;
        n=0;
        tol=0;
        accepted=0;
        step_size=1;
        }
    
        int flag=0;
        while (flag<2)
        {
            x=myran.doub();
            y=-log(x);
            z=myran.doub();
            if (z<= exp(  -0.5*pow(y-1,2)    )  )
            {
                z=myran.doub();
                if (z<0.5)
                {
                y=-y;
                }

                if (flag==0)
                {
                  l_trial=fmod(l+y*step_size,360);
                  if (l_trial<0)
                  {
                  l_trial=0;
                  }
                  
                  
                }
                else
                {
                b_trial=(b+y*step_size);
                if (b_trial>=90)
                {
                 b_trial=90;
                }
                if (b_trial<=-90)
                {
                 b_trial=-90;
                }
                }
                flag=flag+1;
            }
        }
    // evaluate f for this new point
    Hvar hvar(m_data,m_v_ref,m_l_ref,m_b_ref);
    hvar.set_direction(v,l_trial,b_trial);
    hvar.calc_hvar(minr,nb);
  
    if (use_chi2_a==1)
    {
    f_new=hvar.m_chi2_a;
    }
    else
    {
    f_new=hvar.m_chi2_b;
    }

    // accept or reject new point
    accept=0;
    if (f_new<=f)
    {
      accept=1;
    }
    else if (myran.doub()<tol)
    {
      accept=1;
    }
    
    if (accept==1)
    {
      f=f_new;
      l=l_trial;
      b=b_trial;
      //cout << "fp = " << fp << " (l,b) = " << l << " " << b << endl;
      accepted=accepted+1;
      //accept_ratio=float(accepted)/float(n);
      
     // cout<< "\r" << "Acceptence ratio = " << accept_ratio;
      //std::cout << std::flush;
      
      
      if (burn_in==1)
        {
        b_sum=b+b_sum;
        l_sum=l+l_sum;
        k=k+1;
        }
    }
    
    if (n==2000)
    {
    
    m_l_conv=l_sum/k;
    m_b_conv=b_sum/k;
    converged=1;
    //cout<< "\n Scanner converged at (l,b) = " << m_l_conv << " " << m_b_conv << endl;
    //cout << "f = " << f <<endl;
    m_min_chi_a=hvar.m_chi2_a;
    m_min_chi_b=hvar.m_chi2_b;
    }
}

}





void Scanner::rej_method()
{
    ofstream numb;
    numb.open ("numbers.dat");
    int seed;
    double x,y,z;
    seed=5;
    Ran myran(seed);
    for ( int a=1; a<1000; a=a+1){
        x=myran.doub();
        y=-log(x);
        z=myran.doub();
        if (z<= exp(  -0.5*pow(y-1,2)    )  )
        {
        
            z=myran.doub();
            if (z<0.5)
            {
            y=-y;
            }
        
            numb << y <<endl;
        }
    }
    numb.close();
}

