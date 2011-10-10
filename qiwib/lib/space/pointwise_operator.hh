#ifndef POINTWISE_OPERATOR_H
# define POINTWISE_OPERATOR_H

#define inplace_operator(op) \
gridfunction_member(gridfunction<space>&) operator op##=(const gridfunction<space>& g)\
{ for(unsigned int i=0;i<this->size();i++) (*this)[i] op##= g[i]; return *this; }

#define const_operator(op) \
gridfunction_member(gridfunction<space>) operator op (const gridfunction& g) const \
{ gridfunction z(*this); return (z op##= g); }

#define scalar_inplace_operator(op) \
gridfunction_member(gridfunction<space>&) operator op##=(const gridfunction::scalar_t& s) \
{  for(unsigned int i=0;i<this->size();i++) (*this)[i] op##= s; return *this; } \

#define scalar_const_operator(op) \
gridfunction_member(gridfunction<space>) operator op (const scalar_t& s) const \
{ gridfunction z(*this); return (z op##= s); }

#define pointwise_operator(op)\
  inplace_operator(op) \
  const_operator(op)

#define scalar_operator(op) \
  scalar_inplace_operator(op) \
  scalar_const_operator(op)

#endif
