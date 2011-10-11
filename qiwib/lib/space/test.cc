#include <cmath>
#include <complex> 

#include "space.h"
#include "basisset.hh"

using namespace std;

int main()
{
  const int Nx = 1000000;

  RealGrid1D    Rline(-M_PI,M_PI,Nx);
  ComplexGrid1D Cline(-M_PI,M_PI,Nx);

  gridfunction<RealGrid1D>    f(Rline);
  gridfunction<ComplexGrid1D> g(Cline);

  for(int i=0;i<Nx;i++){
    f[i] = sin(2*i*M_PI/((double)Nx));
    //    g[i] = std::exp(complex<double>(0,2)*complex<double>(i*M_PI/((double)Nx),0));
    g[i] = std::exp(complex<double>(0,2*i*M_PI/((double)Nx)));
  }

  gridfunction<RealGrid1D> df(Rline.derivative(f,true));
  gridfunction<ComplexGrid1D> dg(Cline.derivative(g));

  complex<double> int_g  = Cline.integrate(g);
  complex<double> int_dg = Cline.integrate(dg);

  printf("int(f(x)^2)  = % .7g\n",Rline.integrate(f*f));
  printf("int(f(x))    = % .1g\n",Rline.integrate(f));
  printf("int(f'(x)^2) = % .7g\n",Rline.integrate(df*df));
  printf("int(g(x))    = % .1g + % .1g * i\n",int_g.real(), int_g.imag());
  printf("int(g'(x))   = % .1g + % .1g * i\n",int_dg.real(), int_dg.imag());

  basisset<gridfunction<RealGrid1D> > basisset(Rline,10);
  
  return 0;
}
