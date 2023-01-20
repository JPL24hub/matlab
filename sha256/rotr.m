
function [y] = rotr(x, a)
%rotate right

    y = zeros(1,length(x));
    y(1,1:a) = x(1,(length(x)-a+1):end);
    
    for i=1 : length(x)-a
        y(1,a+1) = x(1,i);
        a = a+1;
    end
end