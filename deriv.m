function dydx = deriv(u, k)
% DERIV calculates the first derivative in Fourier space

  dydx = ifft(1i*k.*fft(u));
      
end
