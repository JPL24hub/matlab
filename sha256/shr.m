function [y] = shr(x, a)
%right shift

    y = zeros(1,length(x));
    
    for i=1 : length(x)-a
        y(1,a+1) = x(1,i);
        a = a+1;
    end
end