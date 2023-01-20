function[a] = choice(x, y, z) % if x is 1 take y if x is 0 take z

    a = zeros(1, length(x));
    
    for i=1 : length(x)
        if(x(1,i) == 1)
            a(1,i) = y(1,i);
        else
            a(1,i) = z(1,i);
        end
    end

end