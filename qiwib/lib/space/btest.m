#!/usr/bin/octave
space -global

space = realgrid(-pi,pi,1000);
xs    = space.get_xs();

phi   = realbasis(space,10);

for i = [1:10]
  phi(i-1) = realfunction(sin(i*xs));
end

phi.overlap_matrix()
