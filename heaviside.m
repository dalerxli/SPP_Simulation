function [ y ] = heaviside( x )

%HEAVISIDE is the heaviside function evaluated at x

if x < 0
    y =0;
elseif x > 0
        y = 1;
elseif x == 0
        y=1/2;
end
    
end

