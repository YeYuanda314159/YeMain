function y = yfun(k1,k2,p,x)
%Gong Ye family 
y = ((k1^p-k2^p)*x+k2^p).^(1/p);
end