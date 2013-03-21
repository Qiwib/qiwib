#!/usr/bin/octave
space -global

function newbasis = symmetric_orthogonalize(basis)
  S = basis.overlap_matrix();
  C = S**(-1/2);

  newbasis = basis*C;		# Is actually left multiplication :-/
endfunction

space = realgrid(0,20,10000);
xs    = space.get_xs();

Nb = 5;
phi    = realbasis(space,Nb);

for i = [1:Nb]
  f = realfunction((xs.**i) .* exp(-xs));
  f *= 1/sqrt(space.inner(f,f));
  phi(i-1) = f;
end

printf "Overlap matrix:\n";
phi.overlap_matrix()

printf "Plotting functions:\n";
#subplot(1,2,1);
plot (xs, phi(0).get_data(),xs, phi(1).get_data(),xs, phi(2).get_data(),xs, phi(3).get_data(),xs, phi(4).get_data())

pause

phi = symmetric_orthogonalize(phi);

printf "Orthogonalized basis overlap matrix:\n";
phi.overlap_matrix()

printf "Plotting orthogonalized functions:\n";
#subplot(1,2,2);
plot (xs, phi(0).get_data(),xs, phi(1).get_data(),xs, phi(2).get_data(),xs, phi(3).get_data(),xs, phi(4).get_data())
pause
input('Press Enter to quit!');
