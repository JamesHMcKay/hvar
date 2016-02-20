#ifndef GET_DATA_H
#define GET_DATA_H
#include <vector>
using namespace std;

class Get_data
{
  public:
  std::vector<double> CMBcz;
  std::vector<double> d;
  std::vector<double> CMBvpec;
  std::vector<double> sig;
  std::vector<double> gl;
  std::vector<double> gb;
  int np;
  
  Get_data() {} //constructor
  
  void get_COMPOSITE();
  
  void get_CF2();
  
};

#endif