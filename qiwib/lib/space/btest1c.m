#!/usr/bin/env octave
space -global

space = complexgrid(-pi,pi,2001);
xs    = space.get_xs();
xmin  = space.xmin 

phi   = complexbasis(space,4);

for i = [1:4]
  phi(i-1) = complexfunction(exp(I*i*xs)/sqrt(2*pi));
end

plot(xs,real(phi(0).get_data()),xs,imag(phi(0).get_data()))
pause

S = phi.overlap_matrix()		# Should be identity matrix
T = -0.5*phi.laplacian_matrix()         # Should be -1/2 (j^2 * \delta_{ij})

dphi = phi.derivative();

plot (xs,real(dphi(0).get_data()),xs,imag(dphi(0).get_data()))
pause;


