function[y] = x_or(a, b)

    y = zeros(1, length(a));

    for i=1 : length(a)
        if (a(1,i) == b(1,i))
            y(1,i) = 0;
        else
            y(1,i) = 1;
        end 
    end
end