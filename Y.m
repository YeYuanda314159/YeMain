function y = Y(k1,k2,chi,p)
if p == 0
    y = k1^(chi)*k2^(1-chi);
else
    y = ((k1^p-k2^p)*chi+k2^p).^(1/p);
end