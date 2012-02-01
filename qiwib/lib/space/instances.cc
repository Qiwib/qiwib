#ifndef INSTANCES_CC
# define INSTANCES_CC
#include <complex>
#include "space.h"

#include "grid1d.cc"
#include "basisset.hh"
using namespace std;

template class Grid1D<double>;
template class Grid1D< std::complex<double> >;
template class gridfunction< Grid1D<double> >;
template class gridfunction< Grid1D< complex<double> > >;
//typedef gridfunction< Grid1D<double> > realfunction;
//typedef gridfunction< Grid1D< complex<double> > > complexfunction;

template class basisset<realfunction>;
template class basisset<complexfunction>;

template class std::vector<double>;
template class std::vector< std::complex<double> >;

#endif



