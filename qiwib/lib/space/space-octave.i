#pragma SWIG nowarn=509,124
%include <typemaps.i>
%include <std_complex.i>
%include <stl.i>

%template(vectord) std::vector<double>;

%typemap(in) (const double *v, unsigned int n) 
{
    const octave_value mat_feat=$input;
    if (!mat_feat.is_matrix_type() || !mat_feat.is_double_type() || mat_feat.rows()!=1)
      SWIG_fail;

    const Array<double>& m(mat_feat.vector_value());
    $1 = (double*)m.fortran_vec();
    $2 = mat_feat.columns();
}

%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER)
  (const double *v, unsigned int n)
{
  const octave_value& m($input);

  $1 = (m.is_matrix_type() && m.is_double_type() && m.rows()==1) ? 1 : 0;
}


%typemap(out) std::vector<double> 
{
  unsigned int n = $1.size();
  Matrix mat(dim_vector(1, n));
  
  for (unsigned int i=0; i<n; i++)
    mat(i) = $1[i];

  $result=mat;
}


%typemap(in) (const double *v, unsigned int m, unsigned int n) 
{
    const octave_value mat_feat=$input;
    if (!mat_feat.is_matrix_type() || !mat_feat.is_double_type())
      SWIG_fail;

    const Matrix<double>& m(mat_feat.matrix_value());
    $1 = (double*)m.fortran_vec();
    $2 = mat_feat.rows();
    $3 = mat_feat.columns();
}

%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER)
  (const double *v, unsigned int m, unsigned int n)
{
  const octave_value& m($input);

  $1 = (m.is_matrix_type() && m.is_double_type()) ? 1 : 0;
}


%typemap(out) Array2D<double> 
{
  const size_t m = $1.rows(), n = $1.columns();
  Matrix mat(dim_vector(2,m,n));
  
  for (size_t i=0; i<m; i++)
    for (size_t j=0; j<n; j++)
      mat(i,j) = $1(i,j);

  $result=mat;
}


%include "space.i"

