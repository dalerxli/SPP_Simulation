function dydx = inverse_deriv(u, k)

% INVERSE_DERIV calculates d_x^{-1} in Fourier space

k = [0,  1./(1i*k(2:length(k))) ];

  dydx = ifft(k.*fft(u));
  
%   
% tmp = fft(u)./k;
% tmp(1) = 0;
% 
% dydx = ifft(tmp);  johns mistake?
  
end
