%include "space.i"

%{
#include <octave/config.h>
#include <octave/ov.h>
#include <octave/octave.h>
#include <octave/oct-types.h>
#include <octave/oct-obj.h>
#include <octave/variables.h>
#include <octave/Array.h>
%}

%typemap(typecheck,precedence=SWIG_TYPECHECK_POINTER) (Array)
{
  const octave_value v = $input;
  $1 = (v.is_matrix_type() && v.oct_type_check() && v.rows()==1) ? 1:0;
}

%typemap(in) (const Array&)
{
  fprintf(stderr,"double* input conversion.\n");
  const octave_value v = $input;

  /* if(!v.is_real_matrix()) {  */
  /*   fprintf(stderr,"Input error in passing vector from Octave to SWIG.\n"); */
  /*   SWIG_fail; */
  /* } else { */
  /*   fprintf(stderr,"Input is a row vector.\n"); */
  /* } */

  $1 = new Array<double>(v.vector_value());
};
