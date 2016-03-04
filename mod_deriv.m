function dydx = mod_deriv(u, k)
% MOD_DERIV calculates |D| in Fourier space
  dydx = ifft(abs(k).*fft(u));
end
