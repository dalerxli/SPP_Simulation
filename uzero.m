function u0 = uzero(x,lx, ic_name)
% UZERO defines initial condition  u0 = u(x,0)

amp0 = 1.;

switch ic_name
    case 'sine'
        u0 = amp0*sin(x);
    case 'triangle'
        u0(1:lx/4) = x(1:lx/4)/pi;                           % 0 <= x <= pi/2
        u0((lx/4)+1:3*lx/4) = (pi - x((lx/4)+1:3*lx/4))/pi;  % pi/2 < x <= 3*pi/2
        u0((3*lx/4)+1:lx) = (x((3*lx/4)+1:lx) - 2*pi)/pi;    % 3*pi/2 < x < 2*pi
        u0 = 2.*pi*u0;
     case 'riemann'
        lx0 = round(lx/4);
        lx1 = lx - lx0;
        u0(1:lx0) = 0.;             % 0 <= x <= pi/2
        u0(lx0+1:lx1) = amp0;       % pi/2 < x <= 3*pi/2
        u0(lx1+1:lx) = 0.;          % 3*pi/2 < x < 2*pi
    case 'sawtooth'
        u0(1:lx/2) = x(1:lx/2)/pi;                          % 0 <= x <= pi
        u0((lx/2)+1:lx) = (x((lx/2)+1:lx) - 2*pi)/pi;       % pi < x < 2*pi
    
    
    case 'two_cosines'
        u0 = 2*cos(x ) + 2*cos(2*(x  +2*pi^2));
        
                
     case 'hunter'
         amp1 = 1; 
         amp2 = 2;
         psi0 = amp1.*exp(1i.*x) + amp2.*exp(2i.*( x + 2 * pi^2));
         u0 = psi0 ;%+ conj(psi0);
         
         
   case 'three_cosines'
        u0 = 2*cos(x) + cos(2*(x +2*pi^2)) + cos( 3*(x + 3*pi^3)); 

    case 'fourier_cfs'
        u0 = u0_fft(lx);
    case 'realcf'
        m=2;
        a = zeros(lx/2-1,1);
        for i =1:1:m
        a(i,1) = 20*(rand(1,1)-.5);
        end
        ap = [0, a', 0, fliplr(a')]';
        b = zeros(lx/2-1,1);
        for i = 1:1:m
        b(i,1) = 20*(rand(1,1)-.5);
        end
        bp = [0, b', 0, -fliplr(b')]';
        d = ap + 1i.*bp;
        d = ifft(d);
        u0= d'; 
     case 'imagcf'
        m=4;
        a = zeros(lx/2-1,1);
        b = zeros(lx/2-1,1);
        for i =1:1:m
        a(i,1) = 20*(rand(1,1)-.5);
        b(i,1) = 20*(rand(1,1)-.5);
        end
        ap = [0, a', 0, b']';
        c = zeros(lx/2-1,1);
        d = zeros(lx/2-1,1);
        for i =1:1:m
        c(i,1) = 20*(rand(1,1)-.5);
        d(i,1) = 20*(rand(1,1)-.5);
        end
        %cp = [0, c', 0, fliplr(c')]';
        cp = [0, c', 0, d']';
        d = ap + 1i.*cp;
        d = ifft(d);
        u0= d';
      case 'gaussian'
        a = exp( -20*(x-pi).^2); 
        b = fft(a);
        b(1,1) = 0;
        u0 = ifft(b);
      case 'imaginarygaussian'
        a = 1i.*exp( -20*(x-pi).^2); 
        b = fft(a);
        b(1,1) = 0;
        u0 = ifft(b);

      case 'asym'
         u0 = exp(x).*(2*pi -x).*x;
        
      case 'randcf'
                
        m=10;
        a = zeros(lx/2-1,1);
        ab = zeros(lx/2-1,1);
        for i =1:1:m
        a(i,1) = 20*rand-10;
        end
        for i =1:1:m/2
        %ab(i,1) = rand-.5;
        end
        ap = [0, a', 0, fliplr(ab')]';
        ap = ifft(ap);
        u0= ap'; 
        

         

      case 'randlen'
                
        m=10;
        a = zeros(lx/2-1,1);
        ab = zeros(lx/2-1,1);
        b = zeros(lx/2-1,1);
        bb = zeros(lx/2-1,1);
        
        one = floor(50*rand);
        two = floor(50*rand);
        
        top = max(one,two)+1;
        bottom = min(one,two)+1;
        
        for i =bottom:1:top
        a(i,1) = rand-.5;
        end
        
        one = floor(50*rand);
        two = floor(50*rand);
        
        top = max(one,two)+1;
        bottom = min(one,two)+1;
        
        for i =bottom:1:top
        ab(i,1) = rand-.5;
        end
        
        one = floor(50*rand);
        two = floor(50*rand);
        
        top = max(one,two)+1 ;
        bottom = min(one,two)+1;
        
        for i =bottom:1:top
        b(i,1) = rand-.5;
        end
        
        one = floor(50*rand);
        two = floor(50*rand);
        
        top = max(one,two)+1 ;
        bottom = min(one,two)+1;
        
        for i =bottom:1:top
        bb(i,1) = rand-.5;
        end
        
        ap = [0, a', 0, fliplr(ab')]';
        bp = [0, b', 0, fliplr(bb')]';
        ap = ifft(ap + 1i* bp);
        u0= ap'; 

    case 'poscf'

        m=5;
        a = zeros(lx/2-1,1);
        ab = zeros(lx/2-1,1);
        for i =1:1:m
        a(i,1) = 20*rand-10;
        end
%         for i =1:1:m/2
%         %ab(i,1) = rand-.5;
%         end
        ap = [0, a', 0, fliplr(ab')]';
        ap = ifft(ap);
        u0= ap'; 
    
        
    case 'negcf'
    rng(1)
          m=5;
        a = zeros(lx/2-1,1);
        ab = zeros(lx/2-1,1);
        for i =1:1:m
        a(i,1) = 20*rand-10;
        end
%         for i =1:1:m/2
%         %ab(i,1) = rand-.5;
%         end
        ap = [0, ab', 0, fliplr(a')]';
        ap = ifft(ap);
        u0= ap'; 
    
        
    case 'compact'        
        a = exp( -20*((x-pi)).^2); 
        
        for i= 1:1:lx/4
           a(1,i) = 0;
           a(1,lx-i+1)=0;
        end
        u0 = a;
        
        
        
            case 'compact'        
        a = exp( -20*((x-pi)).^2); 
        
        for i= 1:1:5*lx/11
           a(1,i) = 0;
           a(1,lx-i+1)=0;
        end
        u0 = a;
        
        
                
        
            case 'compact2'        
        a = exp( -20*((x-pi)).^2).*sin(4*x).*cos(-2*x); 
        
        for i= 1:1:lx/4
           a(1,i) = 0;
           a(1,lx-i+1)=0;
        end
        u0 = a;
        
        
         
            case 'compact3'        
        a = 1-5*(x-pi).^2; 
        for i = 1:1:length(a)
            if a(i) < 0
                a(i) = 0;
            end
        end
        
        
        u0 = a;
        
        
        
    otherwise
        disp('No IC specified')
end

u0 = u0.';
end

function u0=u0_fft(lx)
%U0_FFT
    af = zeros(1,(lx/2-1));
    for ii = 1:(lx/2-1)
        af(ii)= cf(ii,lx);
    end
    bf = fliplr(conj(af));
    u0_hat = lx*[0 af 0 bf];
    u0 = ifft(u0_hat);
end

function cn=cf(ii,lx)
   %r1 = rand/sqrt(ii);
   r1 = 1/sqrt(ii); 
   r2 = 2*pi*rand;
   cn = r1*exp(1i*r2);
end