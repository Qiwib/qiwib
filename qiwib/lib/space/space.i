%module space
%include <std_complex.i> 
%include <stl.i> 

%{
#include <complex>
#include "space.h"
%}

%template(vectord)       std::vector<double>;
%template(vectorc)       std::vector< std::complex<double> >;

%include "grid1d/grid1d.hh"

%template(realgrid)        Grid1D<double>;
%template(complexgrid)     Grid1D< std::complex<double> >;

%include "grid1d/grid1d.hh"

%template(realfunction)    gridfunction< Grid1D<double> >;
%template(complexfunction) gridfunction< Grid1D<std::complex<double> > >;


