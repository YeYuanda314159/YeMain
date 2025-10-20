function dy = dyfun(k1,k2,p,x)
%Derivative of Gong Ye family 
dy = 1/p*(k1^p-k2^p)*((k1^p-k2^p)*x+k2^p).^(1/p-1);
end