function dydx = inverse_mod(u, k)

% INVERSE_MOD calculates d_x^{-1} in Fourier space

k = [0,  1./(abs(k(2:length(k))))' ]';

  dydx = ifft(k.*fft(u));
end
