function dydx = deriv2(u, k)

% DERIV2 calculates the laplacian on high wave numbers, smoothly transitions,


l = length(k);
m = max(k);

sq = round(4*sqrt(m));

sq1 = round(5*sqrt(m));

%damp out wave numbers smoothly above sq  below sq1 then damps above sq1
%fully

 visc = fft(u);
 



d = sq1 - sq;


%0 000000sq  . . . . sq1 11 1 1 1 1
% 

k(1:sq+1) = 0;
 k(sq+1:sq1+1) = k(sq+1:sq1+1).*(0:1/d:1)';
 
k(l-sq +1 :l) = 0;
 k(l-sq1+1 : l-sq+1) = k(l-sq1+1: l-sq+1).*fliplr((0:1/d:1))';




%2nd order
dydx = ifft(- (abs(k).^2).*visc);

%fourth order
%dydx = ifft( -(abs(k).^4).*visc);
end

