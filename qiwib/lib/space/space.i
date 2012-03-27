%module space
%{
#include <complex>
#include "space.h"
%}
%include <std_vector.i>
%include "grid1d.hh"

%template(realgrid)        Grid1D<double,double>;
%template(complexgrid)     Grid1D< std::complex<double>, std::complex<double> >;

%include "grid1d.hh"

%template(realfunction)    gridfunction< Grid1D<double, double> >;
%template(complexfunction) gridfunction< Grid1D<std::complex<double>, std::complex<double> > >;

%include "spingrid1d.hh"

* %template(realspin2grid)        SpinGrid1D<double,2>;
* %template(complexspin2grid)     SpinGrid1D< std::complex<double>, 2 >;

* %include "spingrid1d.hh"
* 
* %template(realspin2function)    spingridfunction< double,2 >;
* %template(complexspin2function) spingridfunction< std::complex<double>, 2 >;

%include "grid-instances.hh"

%template(vectorrf) std::vector< gridfunction< Grid1D<double,double> > >;
%template(vectorcf) std::vector< gridfunction< Grid1D<std::complex<double>, std::complex<double> > > >;

%include "basisset.hh"

%template(realbasis) basisset< gridfunction< Grid1D<double,double> > >;
%template(complexbasis) basisset< gridfunction< Grid1D<std::complex<double>, std::complex<double> > > >;
* %template(realspin2basis) basisset< spingridfunction< double,2> >;
* %template(complexspin2basis) basisset< spingridfunction< std::complex<double>, 2 > >;
