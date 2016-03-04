
function yt = f_uv(~,y,k,lx,alpha,beta,nu) 
% F_UV is right hand side of evolution equation psi_t = fuv(psi)

y = y.';
u = y(1:lx);
v = y(lx+1:2*lx);

ut = -1i* P(deriv( alpha*inverse_deriv(u.^2,k).*conj(u)...
       + beta*v.*inverse_deriv(u.*conj(v),k),k) ,lx);
%   ut = -1i* P( alpha*u.^2.*conj(u),lx) ;
vt =- 1i* Q(deriv(alpha*inverse_deriv(v.^2,k).*conj(v)...
      + beta*u.*inverse_deriv(v.*conj(u),k),k),lx);
  
ut = ut - 1i*nu*inverse_deriv(u,k);
vt = vt - 1i*nu*inverse_deriv(v,k);
ut = truncate(ut,lx);
vt = truncate(vt,lx);

yt = [ut, vt];
yt = yt.';
end


