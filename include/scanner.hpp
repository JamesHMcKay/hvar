

#ifndef SCANNER_H
#define SCANNER_H
#include <vector>
using namespace std;

typedef unsigned long long int Ullong;
typedef unsigned int Uint;

class Scanner
{
  private:
  
  double min_b;
  double m_v_ref;
  double m_l_ref;
  double m_b_ref;
  int m_use_H_bar;
  Boost_offset minimum_boost;
  Get_data m_data;
  Hvar m_hvar1;
  
  public:
  
  double min_l;
  double m_min_chi_a;
  double m_min_chi_b;
  double m_l_conv;
  double m_b_conv;
  double m_l_sig;
  double m_b_sig;
  
  Scanner(Hvar hvar1)
  {
  m_v_ref=hvar1.Get_v_ref();
  m_l_ref=hvar1.Get_l_ref();
  m_b_ref=hvar1.Get_b_ref();
  m_data=hvar1.m_data;
  m_hvar1=hvar1;
  m_use_H_bar=0;
  }  // defualt constructor
  
  Scanner(Hvar hvar1,int flag)
  {
  m_v_ref=hvar1.Get_v_ref();
  m_l_ref=hvar1.Get_l_ref();
  m_b_ref=hvar1.Get_b_ref();
  m_data=hvar1.m_data;
  m_hvar1=hvar1;
  m_use_H_bar=flag;
  }  // constructor with option for using model indepdent boost magnitude calculation
  

  
  
  Scanner(Get_data data,double v,double l,double b)
  {
  m_data=data;
  m_v_ref=v;
  m_l_ref=l;
  m_b_ref=b;
  }// alternative constructor for when only interested in chi_a or chi_b (no reference hvar object needed)
  
  
  void find_best_fit(double v_init,double l_init,double b_init,int burnin_length, int chain_length);
  
  void find_min_chi(double v_init,double l_init,double b_init,double minr,int nb,int use_chi2_a);
  
  void find_min_chi_a(double v_init,double l_init,double b_init,double minr,int nb)
  {
  find_min_chi(v_init,l_init,b_init,minr,nb,1);
  };
  
  void find_min_chi_b(double v_init,double l_init,double b_init,double minr,int nb)
  {
  find_min_chi(v_init,l_init,b_init,minr,nb,0);
  };
  
  
  
  
  // return information about the best fit in Boost_offset object
  Boost_offset Get_boost_object() { return minimum_boost; };

  void rej_method(); // generate random number from normal distribution

  struct Ran{
  Ullong v;
  Ran(Ullong j) : v(4101842887655102017LL) {
  v ^= j;
  v = int64();
  }
  inline Ullong int64() {
  v ^= v >> 21; v ^= v << 35; v ^= v >> 4;
  return v * 2685821657736338717LL;
  }
  inline double doub() { return 5.42101086242752217E-20 * int64(); }
  inline Uint int32() { return (Uint)int64(); }
  };

};

#endif