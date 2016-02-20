#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>

#include <fstream>

#include "boost_offset.hpp"

using namespace std;



void Boost_offset::calc_boost_offset(Hvar hvar1,Hvar hvar2) //b,p,S,sigma_b,a)
{




std::vector<double> H1;
std::vector<double> H2;
std::vector<double> av_r2_vec;
std::vector<double> sig_av_r2_vec;
std::vector<double> sig_H1;
std::vector<double> sig_H2;
int nb=hvar1.m_nb;




std::vector<double> Del_H(nb); // would be better to define these after calling number of shells from the Hvar object
std::vector<double> sig_Del_H(nb);
std::vector<double> av_r2(nb);
std::vector<double> sig_av_r2(nb);

std::vector<double> mean_r;

mean_r=hvar2.m_av_r;

H1=hvar1.m_H;
H2=hvar2.m_H;
av_r2_vec=hvar2.m_av_r2;
sig_av_r2_vec=hvar2.m_sig_av_r2;

sig_H1=hvar1.m_sig_H;
sig_H2=hvar2.m_sig_H;



//ofstream myfile;
//myfile.open ("data.txt");



if (m_use_H_bar==1)
{
  for (int i=0;i<nb;i++)
  {
  Del_H[i]=(H2[i]-H1[i])*H1[nb-1];
  sig_Del_H[i]=pow( pow(sig_H1[i],2)+pow(sig_H2[i],2) , 0.5)*H1[nb-1]; // maybe should also add sig_H1[nb-1] into here, have tried this though, negligible effect, so I think this is fine
  av_r2[i]=av_r2_vec[i];
  sig_av_r2[i]=sig_av_r2_vec[i];
  }
}
else
{
  for (int i=0;i<nb;i++)
  {
  Del_H[i]=(H2[i]-H1[i]);
  sig_Del_H[i]=pow( pow(sig_H1[i],2)+pow(sig_H2[i],2) , 0.5);
  av_r2[i]=av_r2_vec[i];
  sig_av_r2[i]=sig_av_r2_vec[i];
  }
}



//myfile.close();






int np=nb-1;

int n=0,p=0;



double r_temp,b_1,b_2,sigma_eu;
double X_bar,Y_bar,Z_sum;

double m_sigma_uu,m_sigma_ee,r_1,r_2,r_3,sigma_a,sigma_b,Zx,Zx2,a,S;





if (Del_H[0]<0) // change sign if first point negative beyond one sigma
  {
    for (int i=0; i<np+1;i++)
      {
      Del_H[i]=-Del_H[i];
      }
      p=1;
  }




for (int i=0;i<nb;i++)
{
 if(Del_H[i]>0)
 {
    n=n+1;   // finds n elements that are non-zero
  }
}
//n=n-1;
if (n<2)
{
      n=0;
      if (p==0){p=1;}
      if (p==1){p=0;}
  
      for (int i=0; i<np+1;i++)
      {
       Del_H[i]=-Del_H[i];
      }
      for (int i=0;i<np+1;i++)
      {
         if(Del_H[i]>0)
         {
          n=n+1;  // finds n elements that are non-zero
          }
       }
   // n=n-1;
}


//cout << "n= "<< n<<endl;

std::vector<double> sigma_uu2(n),sigma_ee2(n),x(n),y(n),sigma_uu(n),sigma_ee(n);  // define arrays to carry these n non-zero points
std::vector<double> log_x(n),log_y(n),Z(n),U(n),V(n),r(n),Wy(n),Wx(n),alpha(n);
std::vector<double> x_adjusted(n);
n=0;
for (int i=0;i<nb;i++)
{
 if(Del_H[i]>0)
 {
  y[n]=Del_H[i];
  x[n]=av_r2[i];
  sigma_uu2[n]=sig_av_r2[i];
  sigma_ee2[n]=sig_Del_H[i];
    n=n+1;
   
  }
}


//cout << "n= "<< n<<endl;


for (int i=0;i<n;i++)
{
  sigma_uu[i]=(sigma_uu2[i]/x[i]);
  sigma_ee[i]=(sigma_ee2[i]/y[i]);
}

for (int i=0;i<n;i++)
{
 log_x[i]=log(x[i]);
 log_y[i]=log(y[i]);
}

sigma_eu=0;

for (int i=0;i<n;i++)
{
Wx[i]=1/sigma_uu[i];
Wy[i]=1/sigma_ee[i];
}

double b=-1;

m_sigma_uu=0, m_sigma_ee=0;
r_1=0, r_2=0, r_3=0;

for (int i=0;i<n;i++)
{
m_sigma_uu=sigma_uu[i]+m_sigma_uu;
m_sigma_ee=sigma_ee[i]+m_sigma_ee;
}
m_sigma_uu=m_sigma_uu/n;
m_sigma_ee=m_sigma_ee/n;

for (int i=0;i<n;i++)
{
r_1=(sigma_uu[i]-m_sigma_uu)*(sigma_ee[i]-m_sigma_ee)+r_1;
r_2=pow((sigma_uu[i]-m_sigma_uu),2)+r_2;
r_3=pow((sigma_ee[i]-m_sigma_ee),2)+r_3;
}
r_temp=r_1/(pow(r_2,0.5)*pow(r_3,0.5));

for (int i=0;i<n;i++)
  {
   r[i]=r_temp;
  }


double b_tmp;
int converged=0,k=0;
while (converged==0)
{
  k=k+1;
  for (int i=0;i<n;i++)
    {
     Z[i]=Wx[i]*Wy[i]/(Wy[i]*(b*b)+Wx[i]-2*b*r[i]*sqrt(Wx[i]*Wy[i]));
    }


  Y_bar=0;
  X_bar=0;
  Z_sum=0;

    for (int i=0;i<n;i++)
      {

      Y_bar=Y_bar+Z[i]*log_y[i];
      X_bar=X_bar+Z[i]*log_x[i];
      Z_sum=Z_sum+Z[i];
      }
  Y_bar=Y_bar/Z_sum;
  X_bar=X_bar/Z_sum;
    for (int i=0;i<n;i++)
      {
      U[i]=log_x[i]-X_bar;
       V[i]=log_y[i]-Y_bar;
       alpha[i]=pow((Wx[i]*Wy[i]),0.5);
      }

  b_1=0;
  b_2=0;
    for (int i=0;i<n;i++)
      {
      b_1=Z[i]*Z[i]*V[i]*(U[i]/Wy[i]+b*V[i]/Wx[i]-r[i]*V[i]/alpha[i])+b_1;
      b_2=Z[i]*Z[i]*U[i]*(U[i]/Wy[i]+b*V[i]/Wx[i]-b*r[i]*U[i]/alpha[i])+b_2;
      }
  b_tmp=b;
  b=b_1/b_2;
  //cout<< "b= " << b << endl;
  
  if (abs(b-b_tmp)<1e-6 && k>2)
  {
  converged=1;
  }
  
}




S=0;
    for (int i=0;i<n;i++)
      {
        S=S+Z[i]* pow ( (log_y[i]-b*log_x[i]-(Y_bar-b*X_bar)),2);
      }


a=exp(Y_bar-b*X_bar);
    for (int i=0;i<n;i++)
      {
       x_adjusted[i]=log_x[i]-((sigma_eu-b*sigma_uu[i])*(log_y[i]-(Y_bar-b*X_bar)-b*log_x[i])/(sigma_ee[i]-2*b*sigma_eu+(pow(b,2))*sigma_uu[i]));
      }
  
  
Zx2=0;
Zx=0;
 for (int i=0;i<n;i++)
      {
Zx=Z[i]*x_adjusted[i]+Zx;
Zx2=Z[i]*pow((x_adjusted[i]),2)+Zx2;
      }
double sigma_beta_0=sqrt(Zx2/(Zx2*Z_sum-pow(Zx,2))); // *a
sigma_b=sqrt(Z_sum/(Zx2*Z_sum-pow(Zx,2)));

sigma_a=sqrt(Zx2/(Zx2*Z_sum-pow(Zx,2)))*a;

upper_A=exp(Y_bar-b*X_bar+sigma_beta_0);
lower_A=exp(Y_bar-b*X_bar-sigma_beta_0);



if (n>1)
{
m_chi_S=S/(n-2+1);  // number of points = n+1 (since we include 0)
}
else
{
m_chi_S=99;
}

m_b=b;






m_A=a;
m_sign_A=p; // 1 if sign has changed
m_sigma_b=sigma_b;
m_sigma_A=sigma_a;



//cout<< "b= " << b << endl;


//cout<< "b , sigma_b, p , a = " << b <<" "<< sigma_b<<" "<< p<<" "<<a<<endl;

//write(90,*) b, sigma_b, p,a


}
