#include <cmath>
#include <complex> 

#include "space.h"
#include "spingrid1d.hh"
#include "basisset.hh"

using namespace std;

int main()
{
  const int Nx = 1000000;

  SpinGrid1D<double,2>        space(-M_PI,M_PI,Nx);
  spingridfunction<double,2>  f(space);

  for(int i=0;i<Nx;i++){
    f[i][0] = cos(2*i*M_PI/((double)Nx));
    f[i][1] = sin(2*i*M_PI/((double)Nx));
  }

  spingridfunction<double,2>  df(space.derivative(f));

  SpinValue<double,2> int_f  = space.integrate(f);
  SpinValue<double,2> int_df = space.integrate(df);
  
  cout << int_f << endl;
  cout << int_df << endl;
  cout << space.inner(f,f) << endl;
  cout << space.inner(df,df) << endl;

  return 0;
}
