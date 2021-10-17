function S_phi = source(x,y)
% Function to calculate the source term
% ft = first term; st = second term
ft = 1000*(2*sinh(x-1/2) + 4*(x-1/2)*cosh(x-1/2)+ (x-1/2)^2*sinh(x-1/2));
st = 1000*(2*sinh(y-1/2) + 4*(y-1/2)*cosh(y-1/2) + (y-1/2)^2*sinh(y-1/2));
S_phi = ft+st;
end

