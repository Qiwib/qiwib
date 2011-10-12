#!/usr/bin/env octave
space -global

space = realgrid(-pi,pi,1000);
xs    = space.get_xs();

phi   = realbasis(space,10);

for i = [1:10]
  phi(i-1) = realfunction(sin(i*xs)/sqrt(pi));
end


S = phi.overlap_matrix()		# Should be pi*I
T = -0.5*phi.laplacian_matrix()

dphi = phi.derivative();

plot (xs,dphi(0).get_data())
pause;


