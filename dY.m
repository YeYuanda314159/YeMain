function y = dY(k1,k2,chi,p)
if p == 0
    y = k1^(chi)*k2^(1-chi)*log(k1/k2);
else
    y = 1/p*(k1^p-k2^p)*((k1^p-k2^p)*chi+k2^p).^(1/p-1);
end