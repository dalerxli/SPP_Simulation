function mdydx = half_deriv(psiq, k)

% HALF_DERIV calculates the half derivative in Fourier space

  mdydx = ifft(sqrt(abs(k)).*fft(psiq));
end
