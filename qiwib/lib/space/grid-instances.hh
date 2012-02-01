#include "grid1d.hh"
#include <algorithm>

typedef Grid1D<double> RealGrid1D;
typedef Grid1D<std::complex<double> > ComplexGrid1D;

class realfunction : public gridfunction<Grid1D<double> > {
public:
  typedef gridfunction<RealGrid1D> baseclass;
  typedef RealGrid1D Space;

  realfunction() {}
  realfunction(const baseclass& b) : baseclass(b) {}
  realfunction(const Space& space) : baseclass(space.Nx) {  }
  realfunction(const double *v, unsigned int n) : baseclass(std::vector<double>(v,v+n)) { 
    fprintf(stderr,"realfunction(pointer) : size(%d)\n",int(this->size()));
    // for(unsigned int i=0;i<this->size();i++)
    //   fprintf(stderr,"(%d %g) ",i, (*this)[i]);
    // fprintf(stderr,"\n");
  }
  realfunction& operator=(const realfunction& b){ 
    resize(b.size());
    copy(b.begin(),b.end(),begin());
    return *this;
  }

};

class complexfunction : public gridfunction<Grid1D<std::complex<double> > > {
public:
  typedef gridfunction<ComplexGrid1D> baseclass;
  typedef ComplexGrid1D Space;

  complexfunction() {}
  complexfunction(const baseclass& b) : baseclass(b) {}
  complexfunction(const Space& space) : baseclass(space.Nx) { }
  complexfunction(const std::complex<double> *v, unsigned int n) : baseclass(std::vector< std::complex<double> >(v,v+n)) { 
    fprintf(stderr,"complexfunction(pointer) : size(%d)\n",int(this->size()));
  }
  complexfunction& operator=(const complexfunction& b){ 
    resize(b.size());
    copy(b.begin(),b.end(),begin());    
    return *this;
  }
};
