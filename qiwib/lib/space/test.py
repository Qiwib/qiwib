#!/usr/bin/env python
from space import realgrid, realfunction;
from numpy import sin,cos,sqrt,pi, linspace, array;
from matplotlib.pyplot import *;

Nx = 10000;
space = realgrid(-pi,pi,Nx);
xs = linspace(-pi,pi,Nx);

sinx = realfunction(sin(xs));
cosx = realfunction(cos(xs));

err  = cosx - space.derivative(sinx);
err2 = sinx + space.second_derivative(sinx);
err3 = cosx + space.third_derivative(sinx);

int_sinx2   = space.integrate(sinx*sinx);
inner_sinx2 = space.inner(sinx,sinx);

print "1D grid with %d points:" % Nx;
print "(cos(x) - d/dx sin(x))     L2-error = % .1g" % sqrt(space.inner(err,err));
print "(sin(x) + d^2/dx^2 sin(x)) L2-error = % .1g" % sqrt(space.inner(err2,err2));
print "(cos(x) + d^3/dx^3 sin(x)) L2-error = % .1g" % sqrt(space.inner(err3,err3));
print "\\int sin(x)   dx = % .1g" % space.integrate(sinx);
print "\\int sin(x)^2 dx = %g (method 1) = %g (method 2) ; minus pi = %.1g,%.1g"  % (int_sinx2,inner_sinx2,
                                                                                     pi - int_sinx2, pi-inner_sinx2);
print "\\int sin(x)^3 dx = % .1g" % space.inner(sinx,sinx,sinx);

fig = figure(figsize=(16,12));
adj = subplots_adjust(hspace=0.4,wspace=0.4);

subplot(2,2,1);
p1 = plot (xs,cosx.get_data());
p2 = plot (xs,space.derivative(sinx));
legend([p1,p2], ["cos(x)","$d/dx\ \\sin(x)$"]);

subplot(2,2,2);
p1 = plot (xs,cosx.get_data());
p2 = plot (xs,-array(space.third_derivative(sinx).get_data()));
legend([p1,p2],["cos(x)","$-d^3/dx^3\ \\sin(x)$"]);

subplot(2,2,3);
p1 = plot (xs,err.get_data());
p2 = plot (xs,err2.get_data());
p3 = plot (xs,err3.get_data());
legend([p1,p2,p3],["$\\cos(x) - d/dx\ \\sin(x)$","$\\sin(x) + d^2/dx^2\ \\sin(x)$","$\\cos(x) + d^3/dx^3\ \\sin(x)$"]);
show();
