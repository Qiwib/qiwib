#ifndef INSTANCES_CC
# define INSTANCES_CC
#include <complex>
#include "space.h"

#include "grid1d/grid1d.cc"

template class Grid1D<double>;
template class Grid1D< std::complex<double> >;
template class gridfunction< Grid1D<double> >;
template class gridfunction< Grid1D< std::complex<double> > >;
template class std::vector<double>;
template class std::vector< std::complex<double> >;

#endif



