#pragma SWIG nowarn=509,124
%include <typemaps.i>
%include <std_complex.i>
%include <stl.i>

%template(vectord) std::vector<double>;
%template(vectorc) std::vector< std::complex<double> >;

%typemap(in) (const double *v, unsigned int n) 
{
    const octave_value mat_feat=$input;
    if (!mat_feat.is_matrix_type() || !mat_feat.is_double_type() 
	|| mat_feat.rows()!=1)
      SWIG_fail;

    const Array<double>& m(mat_feat.vector_value());
    $1 = (double*)m.fortran_vec();
    $2 = mat_feat.columns();
}


%typemap(in) (const std::complex<double> *v, unsigned int n) 
{
    const octave_value mat_feat=$input;
    if (!mat_feat.is_matrix_type() || !mat_feat.is_complex_type() 
	|| mat_feat.rows()!=1)
      SWIG_fail;

    fprintf(stderr,"complex matrix in (v,n)\n");
    const Array<std::complex<double> >& m(mat_feat.complex_vector_value());
    $1 = (std::complex<double>*)m.fortran_vec();
    $2 = mat_feat.columns();
}

%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER)
  (const double *v, unsigned int n)
{
  const octave_value& m($input);

  $1 = (m.is_matrix_type() && m.is_double_type() && m.rows()==1) ? 1 : 0;
}

%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER)
  (const std::complex<double> *v, unsigned int n)
{
  const octave_value& m($input);

  $1 = (m.is_matrix_type() && m.is_complex_type() && m.rows()==1) ? 1 : 0;
}


%typemap(out) std::vector<double> 
{
  unsigned int n = $1.size();
  Matrix mat(dim_vector(1, n));
  
  for (unsigned int i=0; i<n; i++)
    mat(i) = $1[i];

  $result=octave_value(mat);
}

%typemap(out) std::vector<std::complex<double> > 
{
  unsigned int n = $1.size();
  ComplexMatrix mat(dim_vector(1, n));
  
  for (unsigned int i=0; i<n; i++)
    mat(i) = $1[i];

  $result=octave_value(mat);
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

%typemap(in) (const std::complex<double> *v, unsigned int m, unsigned int n) 
{
    const octave_value mat_feat=$input;
    if (!mat_feat.is_matrix_type() || !mat_feat.is_complex_type())
      SWIG_fail;

    const Matrix<std::complex<double> >& m(mat_feat.matrix_value());
    $1 = (std::complex<double>*)m.fortran_vec();
    $2 = mat_feat.rows();
    $3 = mat_feat.columns();
}


%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER)
  (const Array2D<double>&)
{
  const octave_value& m($input);

  $1 = (m.is_matrix_type() && m.is_double_type()) ? 1 : 0;
}

%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER)
  (const Array2D<std::complex<double> >&)
{
  const octave_value& m($input);

  $1 = (m.is_matrix_type() && m.is_complex_type()) ? 1 : 0;
}


%typemap(in) (const Array2D<double>&) 
{
    const octave_value mat_feat=$input;
    if (!mat_feat.is_matrix_type() || !mat_feat.is_double_type())
      SWIG_fail;
    const Matrix& m(mat_feat.matrix_value());
    const double *ptr = m.fortran_vec();
    $1 = new Array2D<double>(mat_feat.rows(),mat_feat.columns(),ptr);
}

%typemap(in) (const Array2D<std::complex<double> >&) 
{
    const octave_value mat_feat=$input;
    if (!mat_feat.is_matrix_type() || !mat_feat.is_complex_type())
      SWIG_fail;

    const ComplexMatrix& m(mat_feat.complex_matrix_value());
    const std::complex<double> *ptr = m.fortran_vec();
    $1 = new Array2D< std::complex<double> >(mat_feat.rows(),mat_feat.columns(),ptr);
}



%typemap(out) Array2D<double> 
{
  const size_t m = $1.rows(), n = $1.columns();
  Matrix mat(dim_vector(m,n));

  for (size_t i=0; i<m; i++)
    for (size_t j=0; j<n; j++){
      mat(i,j) = $1(i,j);
    }
  $result=octave_value(mat);
}

%typemap(out) Array2D<std::complex<double> > 
{
  const size_t m = $1.rows(), n = $1.columns();
  ComplexMatrix mat(dim_vector(m,n));

  for (size_t i=0; i<m; i++)
    for (size_t j=0; j<n; j++){
      mat(i,j) = $1(i,j);
    }
  $result=octave_value(mat);
}


%include "space.i"

