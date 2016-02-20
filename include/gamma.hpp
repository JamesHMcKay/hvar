struct Gamma : Gauleg18 {
//Object for incomplete gamma function. Gauleg18 provides coefficients for Gauss-Legendre quadrature.
//incgammabeta.h
static const int ASWITCH=100;
static const double EPS;
static const double FPMIN;
double gln;
double gammp(const double a, const double x) {
//When to switch to quadrature method. See end of struct for initializations.
//Returns the incomplete gamma function P.a;x/.
if (x < 0.0 || a <= 0.0) throw("bad args in gammp");
if (x == 0.0) return 0.0;
else if ((int)a >= ASWITCH) return gammpapprox(a,x,1); //Quadrature.
else if (x < a+1.0) return gser(a,x); //Use the series representation. else return 1.0-gcf(a,x); Use the continued fraction representation.
}
double gammq(const double a, const double x) {
//Returns the incomplete gamma function Q.a; x/   1   P .a; x/.
if (x < 0.0 || a <= 0.0) throw("bad args in gammq");
if (x == 0.0) return 1.0;
else if ((int)a >= ASWITCH) return gammpapprox(a,x,0); //Quadrature.
else if (x < a+1.0) return 1.0-gser(a,x);// Use the series representation.
else return gcf(a,x); //Use the continued fraction representation.
}
double gser(const double a, const double x)
{
//Returns the incomplete gamma function P.a;x/ evaluated by its series representation. Also sets ln .a/ as gln. User should not call directly.
    double sum,del,ap;
    gln=gammln(a);
    ap=a;
    del=sum=1.0/a;
    for (;;) {
}
}
++ap;
del *= x/ap;
sum += del;
if (fabs(del) < fabs(sum)*EPS)
{
    return sum*exp(-x+a*log(x)-gln);
}


  double gcf(const double a, const double x)
  {
  //Returns the incomplete gamma function Q.a;x/ evaluated by its continued fraction rep- resentation. Also sets ln .a/ as gln. User should not call directly.
  int i;
  double an,b,c,d,del,h;
  gln=gammln(a);
  b=x+1.0-a;
  c=1.0/FPMIN;
  d=1.0/b;
  h=d;
      for (i=1;;i++)
      {

          an = -i*(i-a);
          b += 2.0;
          d=an*d+b;
      //Set up for evaluating continued fraction by modified Lentz’s method ( 5.2) with b0 D 0.
      //Iterate to convergence.

          if (fabs(d) < FPMIN) d=FPMIN;
          c=b+an/c;
          if (fabs(c) < FPMIN) c=FPMIN;
          d=1.0/d;
          del=d*c;
          h *= del;
          if (fabs(del-1.0) <= EPS) break;
      }
  return exp(-x+a*log(x)-gln)*h;
  //Put factors in front.

  }


  double gammpapprox(double a, double x, int psig) {
  //Incomplete gamma by quadrature. Returns P.a;x/ or Q.a;x/, when psig is 1 or 0, respectively. User should not call directly.
  int j;
  double xu,t,sum,ans;
  double a1 = a-1.0, lna1 = log(a1), sqrta1 = sqrt(a1);
  gln = gammln(a);
  //Set how far to integrate into the tail:
  if (x > a1) xu = MAX(a1 + 11.5*sqrta1, x + 6.0*sqrta1); else xu = MAX(0.,MIN(a1 - 7.5*sqrta1, x - 5.0*sqrta1)); sum = 0;
  for (j=0;j<ngau;j++) { //Gauss-Legendre.
             t = x + (xu-x)*y[j];
             sum += w[j]*exp(-(t-a1)+a1*(log(t)-lna1));
          }
          ans = sum*(xu-x)*exp(a1*(lna1-1.)-gln);
  return (psig?(ans>0.0? 1.0-ans:-ans):(ans>=0.0? ans:1.0+ans)); }

double invgammp(double p, double a);
//Inverse function on x of P.a;x/. See  6.2.1.

  double gammln(const double xx) {
  //Returns the value lnŒ .xx/  for xx > 0. int j;
  double x,tmp,y,ser;
  static const Doub cof[14]={57.1562356658629235,-59.5979603554754912, 14.1360979747417471,-0.491913816097620199,.339946499848118887e-4, .465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3, -.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3, .844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5}; if (xx <= 0) throw("bad arg in gammln");
  y=x=xx;
  tmp = x+5.24218750000000000; //Rational 671/128.
  tmp = (x+0.5)*log(tmp)-tmp;
  ser = 0.999999999999997092;
  for (j=0;j<14;j++) ser += cof[j]/++y;
  return tmp+log(2.5066282746310005*ser/x);
  }



};
const double Gamma::EPS = numeric_limits<double>::epsilon(); const double Gamma::FPMIN = numeric_limits<double>::min()/EPS;


struct Gauleg18 {
Abscissas and weights for Gauss-Legendre quadrature.
    static const Int ngau = 18;
    static const double y[18];
    static const double w[18];
};
const double Gauleg18::y[18] = {0.0021695375159141994, 0.011413521097787704,0.027972308950302116,0.051727015600492421, 0.082502225484340941, 0.12007019910960293,
0.16415283300752470, 0.21442376986779355, 0.27051082840644336, 0.33199876341447887, 0.39843234186401943, 0.46931971407375483, 0.54413605556657973, 0.62232745288031077,
 0.70331500465597174, 0.78649910768313447, 0.87126389619061517, 0.95698180152629142};
const double Gauleg18::w[18] = {0.0055657196642445571, 0.012915947284065419,0.020181515297735382,0.027298621498568734, 0.034213810770299537,0.040875750923643261,
0.047235083490265582, 0.053244713977759692,0.058860144245324798,0.064039797355015485,
0.068745323835736408,0.072941885005653087,0.076598410645870640,
0.079687828912071670,0.082187266704339706,0.084078218979661945,
0.085346685739338721,0.085983275670394821};



