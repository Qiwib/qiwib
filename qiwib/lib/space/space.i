%module space
%{
#include <complex>
#include "space.h"
%}

%include "grid1d.hh"

%template(realgrid)        Grid1D<double>;
// %template(complexgrid)     Grid1D< std::complex<double> >;

%include "grid1d.hh"

%template(realfunction)    gridfunction< Grid1D<double> >;
// %template(complexfunction) gridfunction< Grid1D<std::complex<double> > >;


