k1 = 10; k2=1;
x = 0:0.001:1;
p = 2;
q = 0.5;
y1 = 1./((1/k1-1/k2)*x.^q+1/k2);
y2 = x*k1+(1-x)*k2;
y3 = x.^p*k1+(1-x.^p)*k2;
plot(x,y1);
hold on
plot(x,y2);
hold on
plot(x,y3);
legend('y1','y2','y3');