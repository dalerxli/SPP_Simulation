function dydx = hilbert(u, k)
% HILBERT calculates the hilbert transform of u spectrally

  dydx = ifft(-1i.*sign(k).*fft(u));
end
