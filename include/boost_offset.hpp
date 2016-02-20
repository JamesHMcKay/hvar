#ifndef BOOST_OFFSET_H
#define BOOST_OFFSET_H

#include "hvar.hpp"

class Boost_offset
{
private:
  double m_b,m_sign_A,m_sigma_b,m_sigma_b_systematic,m_A,m_sigma_A,m_chi_S;
  int m_use_H_bar;
  public:
  
  double upper_A,lower_A,m_v;
  double Get_b() { return m_b; }
  double Get_sig_b() { return m_sigma_b; }
  double Get_A() {
  double A;
  if (m_use_H_bar==0)
  {
  A= m_A;
  }
  else
  {
  A= m_A*(1/float(70));
  // we can either get "model independent" boost velocities, or values of A to plot with
  // but not both, so if making nice plots of the boost offset it's best to use the option m_use_H_bar=0
  }
  return A;}
  
  double Get_sig_A() {return m_sigma_A;}
  
  double Get_sign_A() { return m_sign_A; }
  double Get_chi_S() { return m_chi_S; }
  double Get_sig_b_systematic() { return m_sigma_b_systematic; }
  double Get_fp()
  {
   double fp;
   if (m_sign_A==0)
    {
    fp=abs(m_b+1);
    }
    else
    {
    fp=2;//2-abs(m_b+1);
    }
  return fp;
  }
  
  double Get_v() {
    double v;
    if (m_use_H_bar==0)
    {
      v=pow(2*m_A*70,0.5);// assume a Hubble constant of 70 km/s/Mpc

    }
    else
    {
     v=pow(2*m_A,0.5);
    }
    return v;
  }
  
  
  
  
  
  
  Boost_offset() { m_use_H_bar=0;} //constructor
  
  Boost_offset(int flag)
  {
  m_use_H_bar=flag;
  }   // alternative constructor to use model independant boost offset magnitude calculation
  
  
  void calc_boost_offset(Hvar hvar1, Hvar hvar2);
  
  void calc_boost_offset_systematic(Hvar hvar1, Hvar hvar2,int nb)
  {
  //double rcut=1.125e2; //!nb=11
  //int nb=hvar1.m_nb;
  //!rcut=1.25d2 !nb=8
  //!shell size for most shells
  //double shell_size=rcut/(nb-2);
  const int n=425; // number of divisions
  double b[n],sb[n],minr,sig_b=0,v[n];
  Get_data data=hvar1.m_data;
  //int nb;
  //nb=hvar1.m_nb;

  for (int i=0;i<n;i++)
  {
  minr=float(i+1)*0.01+2;   //float(i)*(shell_size/float(n));
  hvar1.calc_hvar(minr,nb);
  hvar2.calc_hvar(minr,nb);
  
  calc_boost_offset(hvar1, hvar2);
  b[i]=m_b;
  sb[i]=m_sigma_b;
  v[i]=Get_v();
  //cout<< "b= " << m_b << endl;
  }
  double sum_b=0;
  double max=-99;
  double min=0;
  double vel=0;
  for (int i=0;i<n;i++)
  {
   sum_b=sum_b+b[i];
   
   if(b[i]>max)
   {
   max=b[i];
   }
   if(b[i]<min)
   {
   min=b[i];
   }
   }
   
   

  double av_b=(max-min)/2+min;   //sum_b/n;
  
  double tmp=1;
  for (int i=0;i<n;i++)
  {
  if (abs(av_b-b[i])<tmp)
  {
  tmp=abs(av_b-b[i]);
  sig_b=sb[i];
  vel=v[i];
  }
  }
  
  m_sigma_b=sig_b;
  
  m_sigma_b_systematic=abs((max-min)/2);
  
  m_b=av_b;
  m_v=vel;
  
  // need to add similar code to that used for the central value of the uncertainty to
  // obtain a central value for the boost magnitude
  
  
  }
  
 
//  double Del_H[],
//  double sig_Del_H[],
//  double av_r2[],
//  double sig_av_r2[]
  

  
};



#endif