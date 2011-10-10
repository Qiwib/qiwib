%include <typemaps.i>
%include <std_complex.i>
%include <stl.i>

%{
#include <complex>
%}

%template(vectord) std::vector<double>;
%template(vectorc) std::vector< std::complex<double> >;

%include "space.i"



