#!/usr/bin/env python
from space import realgrid, realfunction, complexgrid, complexfunction;
from math  import sin, cos, pi;

Nx = 10000;
space = realgrid(-pi,pi,Nx);

sinx = realfunction([sin(2*i*pi/Nx) for i in range(0,Nx)]);
cosx = realfunction([cos(2*i*pi/Nx) for i in range(0,Nx)]);

err = cosx - space.derivative(sinx);

int_sinx2   = space.integrate(sinx*sinx);
inner_sinx2 = space.inner(sinx,sinx);

print "1D grid with %d points:" % Nx;
print "(cos(x) - d/dx sin(x)) L2-error = %.1g" % space.inner(err,err);
print "\\int sin(x)   dx = % .1g" % space.integrate(sinx);
print "\\int sin(x)^2 dx = %g (method 1) = %g (method 2) ; minus pi = %.1g,%.1g"  % (int_sinx2,inner_sinx2,
                                                                                     pi - int_sinx2, pi-inner_sinx2);
print "\\int sin(x)^3 dx = % .1g" % space.inner(sinx,sinx,sinx);

