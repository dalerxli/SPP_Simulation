function dydx = Q(u, k)

% Q projects u onto its negative modes
 % v = (u + 1i.*hilbert(u,k))/2;
 n = length(u);
 U = fft(u);
 U( 1 : n/2) = 0;
 U(1) = 0;
 v = ifft(U);
  dydx = v;
end

