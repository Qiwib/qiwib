#!/usr/bin/env octave
space -global

space = realspin2grid(-pi,pi,2000);
xs    = space.get_xs();
xmin = space.xmin 

phi   = realspin2basis(space,4);

for i = [1:4]
  phi(i-1) = realspin2function([cos(i*xs)/sqrt(2*pi);sin(i*xs)/sqrt(2*pi)]);
end

# Components functions (with scalar values) can be extracted by calling component(i)
plot(xs,phi(0).component(0).get_data(),xs,phi(0).component(1).get_data());
pause

S = phi.overlap_matrix()		# Should be identity matrix
T = -0.5*phi.laplacian_matrix()         # Should be -1/2 (j^2 * \delta_{ij})

dphi = phi.derivative();

# ...or the function values can be obtained as a scalar Ngrid x Ncomponent matrix.
plot (xs,dphi(0).get_data())
pause;



