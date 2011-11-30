#!/usr/bin/env octave
space -global

space = realgrid(-pi,pi,2001);
xs    = space.get_xs();
xmin = space.xmin 

phi   = realbasis(space,4);
phi2 = phi;


for i = [1:4]
  phi(i-1) = realfunction(sin(i*xs)/sqrt(pi));
end


S = phi.overlap_matrix()		# Should be pi*I
T = -0.5*phi.laplacian_matrix()

dphi = phi.derivative();

plot (xs,dphi(0).get_data())
pause;


