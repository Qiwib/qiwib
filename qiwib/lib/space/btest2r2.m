#!/usr/bin/octave
space -global

function newbasis = symmetric_orthogonalize(basis)
  S = basis.overlap_matrix();
  C = S**(-1/2);

  newbasis = basis*C;
endfunction

space     = realspin2grid(0,20,10000);
basespace = realgrid(0,20,10000);
xs        = space.get_xs();

Nb = 5;
phi    = realspin2basis(space,2*Nb);

for i = [1:Nb]
  f = realfunction((xs.**i) .* exp(-xs));
  f *= 1/sqrt(basespace.inner(f,f));
  phi(i-1)    = realspin2function([f.get_data(); 0*xs]);
  phi(Nb+i-1) = realspin2function([0*xs; f.get_data()]);
end

printf "Overlap matrix:\n";
phi.overlap_matrix()

printf "Plotting functions:\n";
#subplot(1,2,1);
values = zeros(2*Nb,length(xs),2);
for i = [1:2*Nb]
  values(i,:,:) = phi(i-1).get_data();
end


## First component
plot (xs, values(:,:,1))
pause
## second component
plot (xs, values(:,:,2))
pause 

phi = symmetric_orthogonalize(phi);

printf "Orthogonalized basis overlap matrix:\n";
phi.overlap_matrix()

printf "Plotting orthogonalized functions:\n";
#subplot(1,2,2);

for i = [1:2*Nb]
  values(i,:,:) = phi(i-1).get_data();
end
## First component
plot (xs, values(:,:,1))
pause
## second component
plot (xs, values(:,:,2))
pause 


input('Press Enter to quit!');
