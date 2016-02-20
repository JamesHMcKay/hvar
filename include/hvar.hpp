#ifndef HVAR_H
#define HVAR_H
#include <vector>
#include "get_data.hpp"
#include <boost/math/special_functions/gamma.hpp>
using namespace std;

class Hvar
{
private:
  int m_np_input;
  double m_magnitude;
  double m_l_direction;
  double m_b_direction;
  double m_v_ref;
  double m_l_ref;
  double m_b_ref;
  

  public:
  std::vector<double> m_H;
  std::vector<double> m_sig_H;
  std::vector<double> m_av_r;
  std::vector<double> m_av_r2;
  std::vector<double> m_sig_av_r2;
  std::vector<double> m_delta_H;
  std::vector<double> m_sig_delta_H;
  //std::vector<double> prob;
  double prob;
  int m_nb;
  double m_chi2_a;
  double m_chi2_b;
  Get_data m_data;
  int m_nb_prev;
  double m_minr_a;
  double m_minr_b;
  
  
  double Get_v_ref() { return m_v_ref; }
  double Get_l_ref() { return m_l_ref; }
  double Get_b_ref() { return m_b_ref; }
  
  Hvar() {}  // defualt constructor, only used when declaring an hvar
             // that is later assigned to a pre defined hvar object
  
  Hvar(Get_data data) {
      m_v_ref=318.6;
      m_l_ref=106;
      m_b_ref=-6;
      m_data=data;
      m_minr_a=0;
      m_minr_b=0;
      m_nb_prev=0;
  }  // constructor sets LG as the reference frame
  
  Hvar(Get_data data,double v_ref, double l_ref, double b_ref)
  {
  m_v_ref=v_ref;
  m_l_ref=l_ref;
  m_b_ref=b_ref;
  m_data=data;
  m_minr_a=0;
  m_minr_b=0;
  m_nb_prev=0;
  }  // alternative constructor to define a reference frame (must be wrt to sun)
  // for CMB this is (371,264.14,48.26)
  
  void set_direction(double magnitude, double l_direction, double b_direction)
  {
  m_magnitude=magnitude;
  m_l_direction=l_direction;
  m_b_direction=b_direction;
  }
  
  void calc_hvar_add(double minr, int nb, int nb_prev); // use when adding a second
                                                        // shell configuration
  
  void calc_hvar(double minr, int nb)
  {
  calc_hvar_add(minr, nb, 0);
  }
  
//
//  double incgamma (double x, double a)
//  {
//  double sum=0;
//  double term=1.0/a;
//  int n=1;
//    while (term != 0){
//      sum = sum + term;
//      term = term*(x/(a+n));
//      n++;
//      cout<< "here " << endl;
//      }
//  return pow(x,a)*exp(-1*x)*sum;
//  }
  
  
  
  
};

#endif