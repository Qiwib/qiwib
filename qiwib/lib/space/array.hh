#ifndef ARRAY_HH
# define ARRAY_HH
#include <inttypes.h>
#include <string.h>
#include <vector>

template <typename T> class Array2D : public std::vector<T> {
  size_t m,n;
public:
  Array2D(size_t m, size_t n, const T *src) : std::vector<T>(src,src+m*n), m(m), n(n) {}
  Array2D(size_t m=1, size_t n=1) : std::vector<T>(m*n), m(m), n(n) {}

  inline T& operator()(size_t i, size_t j){ return (*this)[i*n+j]; }
  inline const T operator()(size_t i, size_t j) const { return (*this)[i*n+j]; }

  size_t rows() const { return m; }
  size_t columns() const { return n; }

};
#endif
