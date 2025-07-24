function y = maxcos(x,p)
    y = cos(x)^2;
    for i = 2:p
        if cos(i*x)^2 > y
            y = cos(i*x)^2;
        end
    end
end