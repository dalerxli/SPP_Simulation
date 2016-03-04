function u_t = hamiltonian(u) 

% HAMILTONIAN computes the hamiltonian spectrally, computationally ineffecient, better to use spatial form

lx = length(u);
u = u.';
% computes flux as convolution of fft's
n = lx;
u_hat = fft(u);

%conjugate inside or outside?
u_hatc = conj(fft(u));

u_shift = fftshift(u_hat);
u_shiftc = fftshift(u_hatc);

u_pad = [zeros(1,n) u_shift(2:n) zeros(1,n)];
u_padc = [zeros(1,n) u_shiftc(2:n) zeros(1,n)];


%% Compute hamiltonian of previous step

E = 0;
for l = 1:n 
    ll = l - n/2 - 1 ;     %corresponding wave number            the main one, 1
   for k = 1:n   
        kk = k - n/2 -1 ;  %corresponding wave number      the 4
       for j =  1:n 
            jj = j - n/2 - 1 ;   %corresponding wave number           the 3 

      %em wave interaction
    E = E +  lambdah(ll, jj +kk -ll,jj,kk)* u_shift(k) * u_shift(j) * u_padc(jj + kk - ll + 3*n/2)* u_shiftc(l);
       end
   end
end



u_t = E/n^4;
end

%%
function lcf = lambdah(k1,k2,k3,k4)
% LAMBDAH computes three-wave interaction coefficient
mod = abs(k1) + abs(k2) + abs(k3) + abs(k4);
if (mod==0)
    lcf = 0;
elseif ( k1*k2*k3*k4 == 0)
    lcf = 0;
else
    lcf = 1/2*(2*  (k3*k2 + abs(k3*k2))    *    (k4*k1 + abs(k4*k1))    /     (abs(k1) + abs(k2) + abs(k3) + abs(k4))  ...
        + (k3*k4 - abs(k3*k4))    *     (k2*k1 - abs(k2*k1))    /     (abs(k1) + abs(k2) + abs(k3) + abs(k4)) );
end
end