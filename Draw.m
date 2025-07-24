clc
clear
x = 0:0.0000001:pi/2;
y = zeros(length(x),1);
p = 10;
miny = 1;
for i = 1:length(x)
    y(i) = maxcos(x(i),p);
    if y(i) < miny
        miny = y(i);
        minx = x(i);
        minindx = i;
    end
end
plot(x,y)