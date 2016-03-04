function dydx = P(u, k)
% P projects a functon onto its positive modes

 % v = (u + 1i.*hilbert(u,k))/2;
 n = length(u);
 U = fft(u);
 U( n/2+1 : end) = 0;
 U(1) = 0;
 v = ifft(U);
  dydx = v;
end
