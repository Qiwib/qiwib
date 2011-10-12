#!/usr/bin/env octave
space -global;

tic

Nx = 10000;
space = realgrid(-pi,pi,Nx);
xs = [-pi:space.dx:pi];

sinx = realfunction(sin(xs));
cosx = realfunction(cos(xs));
sinx = 2*sinx
err  = cosx - space.derivative(sinx,true);
err2 = sinx + space.second_derivative(sinx,true);
err3 = cosx + space.third_derivative(sinx,true);



int_sinx2   = space.integrate(sinx*sinx);
inner_sinx2 = space.inner(sinx,sinx);

toc

printf( "1D grid with %d points:\n",Nx);
printf("(cos(x) - d/dx sin(x))     L2-error = % .1g\n", sqrt(space.inner(err,err)));
printf("(sin(x) + d^2/dx^2 sin(x)) L2-error = % .1g\n", sqrt(space.inner(err2,err2)));
printf("(cos(x) + d^3/dx^3 sin(x)) L2-error = % .1g\n", sqrt(space.inner(err3,err3)));

printf("\\int sin(x)   dx = % .1g\n", space.integrate(sinx));
printf("\\int sin(x)^2 dx = %g (method 1) = %g (method 2) ; minus pi = %.1g,%.1g\n", 
      int_sinx2,inner_sinx2, pi - int_sinx2, pi-inner_sinx2);
printf("\\int sin(x)^3 dx = % .1g\n",space.inner(sinx,sinx,sinx));


plot (xs,cosx.get_data(),";cos(x);",xs,space.derivative(sinx,true).get_data(),";d/dx sin(x);")
pause;

plot (xs,cosx.get_data(),";cos(x);",xs,-space.third_derivative(sinx,true).get_data(),";-d^3/dx^3 sin(x);")
pause;


plot (xs,err.get_data(),";cos(x) - d/dx sin(x);",
      xs,err2.get_data(),";sin(x) + d^2/dx^2 sin(x);",
      xs,err3.get_data(),";cos(x) + d^3/dx^3 sin(x);")
pause;

