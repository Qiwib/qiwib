#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric;

template <typename Space, typename Function> class basisset : public std::vector<Function> {
public:
  typedef typename Space::value_t  value_t;
  typedef typename Space::scalar_t scalar_t;
  typedef typename std::vector<value_t>::const_iterator const_iterator;
  typedef typename std::vector<value_t>::iterator iterator; 
  
  typedef ublas::symmetric_matrix<scalar_t,ublas::upper> symmetric_matrix;

  const Space&     space;
  symmetrix_matrix overlap;

  bassisset(const Space& space, std::vector<Function> basis) : space(space), 
							       std::vector<Function>(basis.begin(),basis.end()),
							       overlap(size(),size())
  {
    update_overlap_matrix();
  }

  basisset orthogonalize_symmetric()
  {
    
  }

  void update_overlap_matrix()
  {
    unsigned int i = 0, j =0;

    for(const_iterator f(this->begin()); f!=this->end();f++,i++){
      j = i;
      for(const_iterator g(f); g != this->end; g++,j++)
	overlap(i,j) = space.inner(*f,*g);
    }
  }
};
