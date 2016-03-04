function k = make_k(n)
% MAKE_K indexes wavenumber k for fft correctly

k = zeros(1,n);
k(1:(n/2)-1) = -(heaviside(n/2 + 1 - (1:(n/2)-1)).*(-1*((1:(n/2)-1) - 1)) ...
                   + heaviside((1:(n/2)-1) - n/2).*(n - (1:(n/2)-1) + 1));
k((n/2)+2:n) = -(heaviside(n/2 + 1 - ((n/2)+2:n)).*(-1*(((n/2)+2:n) - 1)) ...
                   + heaviside(((n/2)+2:n) - n/2).*(n - ((n/2)+2:n) + 1));
k(n/2) = n/2 - 1;
k((n/2)+1) = -n/2;
end
