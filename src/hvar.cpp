#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>

#include <boost/math/special_functions/gamma.hpp>
#include <fstream>

#include "boost_offset.hpp"
#include "hvar.hpp"
#include "get_data.hpp"

//#include "standalone_SSM.hpp"
using namespace std;

void galactic(double &vtemp_x,double &vtemp_y, double &vtemp_z ,double &vgl_mag,double &vgl_l,double &vgl_b)
{
// computes velocity in galactic polar coordinates
double deg=1.7453292519943295769e-2;
double Pi = 3.1415927;
double vtemp[3],vgl[3];
vtemp[0]=vtemp_x;
vtemp[1]=vtemp_y;
vtemp[2]=vtemp_z;
vgl[0]=pow(  pow(vtemp[0],2)+ pow(vtemp[1],2)+ pow(vtemp[2],2)   ,0.5);
if(vgl[0]!=0)
{
 vgl[2]=asin(vtemp[2]/vgl[0]);
 if (vtemp[0]!=0)
 {
         if ((vtemp[0]<0) && (vtemp[1]<0))
          {
           vgl[1]=-atan(vtemp[1]/vtemp[0]);
          }
          else
          {
           vgl[1]=atan(vtemp[1]/vtemp[0]);
          }
  }
  else
  {
  vgl[1]=(Pi/2);
  }
}
else
{
vgl[2]=0;
vgl[1]=0;
}

vgl[2]=vgl[2]/deg;
vgl[1]=vgl[1]/deg;
if(vtemp[0]<0)
  {
   if(vtemp[1]<0)
     {
     vgl[1]=180-vgl[1];
     }
     else
     {
     vgl[1]=vgl[1]+180;
     }
  }
  else
  {
   if(vtemp[1]<0)
   {
   vgl[1]=vgl[1]+360;
   }
  }
vgl_mag= vgl[0];
vgl_l=vgl[1];
vgl_b=vgl[2];
}



void cartesian(double &vs, double &l,double &b,double &v_x,double &v_y,double &v_z)
{
// computes velocity in Cartesian coordinatesgalactic polar coordinates
v_x=vs*cos(l)*cos(b);
v_y=vs*sin(l)*cos(b);
v_z=vs*sin(b);
}


void Hvar::calc_hvar_add(double minr, int nb,int nb_prev)
{
int i;
std::vector<int >count(nb);
m_nb=nb+nb_prev;
m_nb_prev=nb_prev;
if (nb_prev!=0)
{
m_minr_b=minr;
}
else
{
m_minr_a=minr;
}


const int n_shells=m_nb;
double deg,maxr,bs,H0,S,Delr,Delr3,Delz3,sigsol,lu;
double vCMB,bCMB,lCMB,sbCMB,cbCMB,vREF,bREF,lREF,sgb,cgb,gli,projCMB;
double v,b,l,sb,cb,proj,Delr22,vg[3];
double sumr;
std::vector<double> H(n_shells), si(n_shells),s(n_shells),dmean2(n_shells);
std::vector<double> dH(n_shells),sdH(n_shells),sigr(n_shells);
double rpen,rcut,penr,vrel2;
std::vector<double> dmean(n_shells),szero(n_shells);
double v3[3],v_REF[3],v_CMB[3];

//! 1 degree in radians
      deg=1.7453292519943295769e-2;
//! CMB relative to sun quoted at NED
      vCMB=371;
//! Fixsen et al 1996
//!      vCMB=370.6d0
      lCMB=264.14*deg;
      bCMB=48.26*deg;
//! Local group velocity coordinates Tully et al 2008 (relative to sun)





vREF=m_v_ref;
lREF=m_l_ref*deg;
bREF=m_b_ref*deg;

v=m_magnitude;
l=m_l_direction*deg;
b=m_b_direction*deg;

cartesian(vCMB,lCMB,bCMB,v_CMB[0],v_CMB[1],v_CMB[2]);
cartesian(v,l,b,v3[0],v3[1],v3[2]);
cartesian(vREF,lREF,bREF,v_REF[0],v_REF[1],v_REF[2]);

for (int k=0;k<3;k++)
{
v3[k]=v3[k]+v_REF[k];
}

//cout<< "input cartesian velocity = " << v3[0]<<" "<<v3[1]<<" "<<v3[2]<<endl;

galactic(v3[0],v3[1],v3[2],vg[0],vg[1],vg[2]);

//cout<< "output (v,l,b) = " << vg[0]<<" "<<vg[1]<<" "<<vg[2]<<endl;

v=vg[0];
l=vg[1]*deg;
b=vg[2]*deg;

//! Sine and cosine of CMB, LG frames
      sbCMB=sin(bCMB);
      cbCMB=cos(bCMB);
      //sbLG=sin(bLG);
      //cbLG=cos(bLG);
      sb=sin(b);
      cb=cos(b);

//!Mean Hubble constant in h Mpc
H0=100;
//!Uncertainty in boost of sun wrt LG from Tully 2008
sigsol=(2.01e1)/H0;
//!Minimum radius in /h Mpc
//!minr=dfloat(min)*6.25d0



// Call data from Get_data class

  std::vector<double> CMBcz=m_data.CMBcz;
  std::vector<double> d=m_data.d;
  std::vector<double> CMBvpec=m_data.CMBvpec;
  std::vector<double> sig=m_data.sig;
  std::vector<double> gl=m_data.gl;
  std::vector<double> gb=m_data.gb;
  const int np=m_data.np;
  



std::vector<double> cz(np),sig2(np);
  
for (int i=0;i<np;i++)
{
   //   d[i]=d[i]*0.745;// ! Only required when using CF2 data
      sgb=sin(deg*gb[i]);
      cgb=cos(deg*gb[i]);
      gli=deg*gl[i];
      projCMB=sgb*sbCMB+cgb*cbCMB*cos(gli-lCMB);
      //projLG=sgb*sbLG+cgb*cbLG*cos(gli-lLG);
      proj=sgb*sb+cgb*cb*cos(gli-l);
      //vrel=vLG*projLG-vCMB*projCMB;
      vrel2=v*proj-vCMB*projCMB;
      cz[i]=CMBcz[i]+vrel2;
     // LGcz[i]=CMBcz[i]+vrel;
    //  sig[i]=sig[i]*d[i]*H0;// ! Only required when using CF2 data
//! Convert velocity variance to distance variance
      sig2[i]=sig[i]*sig[i]/H0/H0;
}



  
//!Penultimate unprimed shell inner radius in /h Mpc
//!rcut=1.375d2 !nb=13
rcut=1.125e2; //!nb=11
//!rcut=1.25d2 !nb=8
//!shell size for most shells
bs=rcut/(nb-2);
// cout<< "shell size is " << bs <<endl;
//!bs=(rcut-minr)/dfloat(nb-2) !nb=8, prime shells
//!Last shell inner radius in /h Mpc
penr=1.5625e2;
//inr[nb]=penr;
//!Size extension in penultimate shell radius with boundary at 156.25/h Mpc
rpen=penr-bs*(nb-1)-minr;
//!Maximum radius in /h Mpc
maxr=4.2e2;
i=0;
lu=minr+bs;
std::vector<double> chin(np),csq(nb),nu(nb),i1(nb+1);
    i1[0]=1;
    for (int k=0;k<nb;k++)
    {
    
    // cout<< "lu = " << lu << endl;
     
     count[k]=0;
     S=0;
     Delr=0;
     Delr3=0;
     Delz3=0;
     Delr22=0;
     sumr=0;
     while (d[i]<=minr)
     {
     i=i+1;
     i1[1]=i;
     }
    // istart=i;
      
     while ((d[i]<=lu) && (i<np))
     {
      count[k]=count[k]+1;
      S=S+float(1)/sig2[i];
      
      Delr=Delr+d[i]/sig2[i];
      Delr22=Delr22+(pow(d[i],2))/sig2[i];
      sumr=d[i]+sumr;
  //! Add (cz/sig)^2
      Delz3=Delz3+cz[i]*cz[i]/sig2[i];
  //! Add (cz*d/sig^2)
      Delr3=Delr3+cz[i]*d[i]/sig2[i];
      
  //!  check shell splitting
     i=i+1;
     }
     dmean[k]=Delr/S;
     dmean2[k]=Delr22/S;
     //! Percentage zero point error
      szero[k]=sigsol/dmean[k];
      H[k]=Delz3/Delr3;
      //! Sigma in H units km/s/Mpc
      si[k]=pow(  Delz3 , 1.5)/Delr3/Delr3;
      s[k]=pow  ( (pow(si[k],2)+pow(szero[k],2)*pow(H[k],2)), 0.5);
      sigr[k]=2*pow((Delr22),0.5)/S;
      //! Determine chi square goodness of fit
      nu[k]=(count[k]-1);

      if (k>nb-4)
      {
        lu=lu+rpen;
        //! Increase rpen for final shell
        rpen=maxr-lu-bs;
      }
      lu=lu+bs;
      m_H.push_back(1);
      m_sig_H.push_back(1);
      m_av_r.push_back(1);
      m_av_r2.push_back(1);
      m_sig_av_r2.push_back(1);
      m_H[k+nb_prev]=H[k];
      m_sig_H[k+nb_prev]=s[k];
      m_av_r2[k+nb_prev]=dmean2[k];
      m_av_r[k+nb_prev]=dmean[k];
      m_sig_av_r2[k+nb_prev]=sigr[k];
     }


  for (int k=1;k<nb+1;k++)
  {
  i1[k]=i1[k-1]+count[k-1];
  }

  for (int k=0;k<nb;k++)
  {
     dH[k]=H[k]/H[nb-1];
     sdH[k]=dH[k]*pow(  ( pow((s[k]/H[k]),2) + pow((s[nb-1]/H[nb-1]),2)) , 0.5);
     dH[k]=dH[k]-1;
     //csq[k]=pow((dH[k]/sd[k]),2);
     m_delta_H.push_back(1);
     m_sig_delta_H.push_back(1);
     m_delta_H[k+nb_prev]=dH[k];
     m_sig_delta_H[k+nb_prev]=sdH[k];
     csq[k]=pow( (dH[k]/sdH[k]) , 2 );

  }
  int j=0;
  double chisq1;
  for (int k=0;k<nb;k++)
  {
  chisq1=0;
  j=i1[k]-1;
      while (j<i1[k+1]-1)
      {
          chisq1=chisq1+pow(   (d[j]-cz[j]/H[k]) , 2  )/sig2[j];
        // cout<< sig2[j]<< endl;
        // cout<< d[j]<< endl;
          j=j+1;
      }
  chin[k]=chisq1;
  }
  
  
  double chi2_a=0;
  for (int k=0;k<8;k++)
  {
  chi2_a=csq[k]+chi2_a;
  }
  double chi2_b=0,nu_total=0;
  for (int k=6;k<nb;k++)
  {
  chi2_b=chin[k]+chi2_b;
//  cout<< nu[k] << endl;
  nu_total=nu[k]+nu_total;
  }
  m_chi2_a=chi2_a;
  m_chi2_b=chi2_b/(nu_total+4);
  
  
  
  double dof=0,chi2=0,chisq=0,p1;
  std::vector<double> inr(nb);
  int kk=0;
  for (int k=nb-2-3; k>-1;k--)
  {
    dof=dof+0.5;

    chisq=chisq+csq[k];
    chi2=chisq/2;
    p1=boost::math::gamma_q(dof,chi2)/ (std::tgamma(dof));
    inr[k]=minr+(k-1)*bs;
    //cout << "p1 = " << p1 << endl;
    prob=p1;
    //prob.push_back(1);
    //prob[kk+nb_prev]=p1;
    kk=kk+1;
  
  }
  
  
  
  
  
  
  
  
}





