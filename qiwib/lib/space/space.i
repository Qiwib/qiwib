%module space
%{
#include <complex>
#include "space.h"
%}
%include <std_vector.i>
%include "grid1d.hh"

%template(realgrid)        Grid1D<double>;
// %template(complexgrid)     Grid1D< std::complex<double> >;

%include "grid1d.hh"

%template(realfunction)    gridfunction< Grid1D<double> >;
// %template(complexfunction) gridfunction< Grid1D<std::complex<double> > >;

%template(vectorrf) std::vector< gridfunction< Grid1D< double > > >;

%include "basisset.hh"

%template(realbasis) basisset< gridfunction< Grid1D<double> > >;
