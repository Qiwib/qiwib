%module space
%{
#include <complex>
#include "space.h"
%}
typedef std::complex<double> complex_t;

%include <std_vector.i>
%include "grid1d.hh"

%template(realgrid)        Grid1D<double,double>;
%template(complexgrid)     Grid1D< std::complex<double> , std::complex<double>  >;

%template(realfunction_)    gridfunction< Grid1D<double, double> >;
%template(complexfunction_) gridfunction< Grid1D<std::complex<double> , std::complex<double>  > >;

%include "spingrid1d.hh"

%template(realspinvalue2_)       SpinValue<double,2>;
%template(rsvector_)             std::vector< SpinGrid1D<double,2>::value_t >;
%template(csvector_)             std::vector< SpinGrid1D<std::complex<double>,2>::value_t >;
%template(realspin2function_)    gridfunction< SpinGrid1D<double,2> >; 
%template(complexspin2function_) gridfunction< SpinGrid1D<std::complex<double>,2> >; 
%template(realspingrid2_)        Grid1D<double,SpinValue<double,2> >;
%template(complexspingrid2_)     Grid1D< std::complex<double> , SpinValue<std::complex<double> ,2> >;

%template(realspin2grid)        SpinGrid1D<double,2>; 
%template(complexspin2grid)     SpinGrid1D< std::complex<double> , 2 >; 
%template(realspin2function)    spingridfunction< double,2 >; 
%template(complexspin2function) spingridfunction< std::complex<double> , 2 >; 


%include "grid-instances.hh"

%template(vectorrf) std::vector< gridfunction< Grid1D<double,double> > >;
%template(vectorcf) std::vector< gridfunction< Grid1D<std::complex<double> , std::complex<double>  > > >;
%template(vectorrf2) std::vector< spingridfunction<double,2> >;
%template(vectorcf2) std::vector< spingridfunction< std::complex<double>, 2 > >;

%include "basisset.hh"

%template(realbasis) basisset< gridfunction< Grid1D<double,double> > >;
%template(complexbasis) basisset< gridfunction< Grid1D<std::complex<double> , std::complex<double>  > > >;
%template(realspin2basis) basisset< spingridfunction< double,2> >; 
%template(complexspin2basis) basisset< spingridfunction< std::complex<double> , 2 > >; 
