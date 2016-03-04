function out = truncate (in,n)
% TRUNCATE truncate n/2 high wave numbers
        in = fft(in); 
        in((n/4+2):3*n/4) = 0;

       % out = real(ifft(in));
        out = ifft(in);


%out  = in;
end
