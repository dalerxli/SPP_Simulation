function k = negative(u)

%NEGATIVE gives u(-x) from u(x)

lx = length(u);
a = flipud(u);
k = [ a(lx)', a(1:lx-1)']';

end
