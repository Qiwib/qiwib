#ifndef INSTANCES_CC
# define INSTANCES_CC
#include <complex>
#include "space.h"

#include "grid1d.cc"
#include "basisset.hh"
using namespace std;

typedef complex<double> complex_t;
template class Grid1D<double>;
template class Grid1D< complex_t >;
template class Grid1D<double, SpinValue<double,2> >;
template class Grid1D<complex_t, SpinValue<complex_t,2> >;
template class gridfunction< Grid1D<double> >;

template class gridfunction< Grid1D< complex_t > >;
template class gridfunction< Grid1D<double, SpinValue<double,2> > >;
template class gridfunction< Grid1D<complex_t, SpinValue<complex_t,2> > >;
//typedef gridfunction< Grid1D<double> > realfunction;
//typedef gridfunction< Grid1D<complex_t> > complexfunction;
//typedef gridfunction< Grid1D<double, SpinValue<double,2> > > realspinfunction;
//typedef gridfunction< Grid1D<complex_t, SpinValue<complex_t,2> > > complexspinfunction;
template class spingridfunction<double,2>;
template class spingridfunction<complex_t,2>;


template class basisset<realfunction>;
template class basisset<complexfunction>;
template class basisset<spingridfunction<double,2> >;
template class basisset<spingridfunction<complex_t,2> >;

template class std::vector<double>;
template class std::vector< complex_t >;

#endif





