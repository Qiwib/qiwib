#ifndef ARRAY_HH
# define ARRAY_HH
#include <string.h>

template <typename T> class Array2D {
  unsigned int m,n;
  T *data;
public:
  
  // Array2D(const Array2D& v) : m(m),n(n),data(new T[m*n]) {}
  // {
  //   if(v.data) memcpy(data,v.data,m*n*sizeof(T));
  // }

  Array2D(size_t m=1, size_t n=1, const T *src = 0) : m(m), n(n), data(new T[m*n]) 
  {
    if(src) memcpy(data,src,m*n*sizeof(T));
  }

  ~Array2D(){ delete data; }

  Array2D get_data() const {
    return * new Array2D(*this);
  }

  size_t rows()    const { return m; }
  size_t columns() const { return n; }

  inline T& operator[](size_t i){ return data[i]; }
  inline const T operator[](size_t i) const { return data[i]; }
  inline T& operator()(size_t i){ return data[i]; }
  inline const T operator()(size_t i) const { return data[i]; }

  inline T& operator()(size_t i, size_t j){ return data[i*m+j]; }
  inline const T operator()(size_t i, size_t j) const { return data[i*m+j]; }
};
#endif
