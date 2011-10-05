#!/usr/bin/env octave
space -global;

Nx = 10000;
space = realgrid(-pi,pi,Nx);
sinx = realfunction(space);
cosx = realfunction(space);

tic
for i = [0:Nx-1]
  sinx(i) = sin(2*pi*i/Nx);
  cosx(i) = cos(2*pi*i/Nx);
end
toc

err = cosx - space.derivative(sinx);

int_sinx2   = space.integrate(sinx*sinx);
inner_sinx2 = space.inner(sinx,sinx);

printf( "1D grid with %d points:\n",Nx);
printf("(cos(x) - d/dx sin(x)) L2-error = % .1g\n", space.inner(err,err));
printf("\\int sin(x)   dx = % .1g\n", space.integrate(sinx));
printf("\\int sin(x)^2 dx = %g (method 1) = %g (method 2) ; minus pi = %.1g,%.1g\n", 
      int_sinx2,inner_sinx2, pi - int_sinx2, pi-inner_sinx2);
printf("\\int sin(x)^3 dx = % .1g\n",space.inner(sinx,sinx,sinx));


