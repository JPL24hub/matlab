function [R] = binary_inversion_algorithm(p, a)
%https://www.allaboutcircuits.com/technical-articles/how-to-vhdl-description-of-a-simple-algorithm-the-data-path/
%https://www.allaboutcircuits.com/technical-articles/how-to-vhdl-description-of-a-simple-algorithm-the-control-path/
    u = a;
    v = p;
    x = 1;
    y = 0;

    while u ~= 0
        while mod(vpi(u), vpi(2)) == 0
            u = vpi(u) / vpi(2);
            if mod(vpi(x), vpi(2)) == 0
                x = vpi(x) / vpi(2);
            else
                x = (vpi(x) + vpi(p)) / vpi(2);
            end
        end
        while mod(vpi(v), vpi(2)) == 0
            v = vpi(v) / vpi(2);
            if mod(vpi(y), vpi(2)) == 0
                y = vpi(y) / vpi(2);
            else
                y = (vpi(y) + vpi(p)) / vpi(2);
            end
        end
        if vpi(u) >= vpi(v)
            u = vpi(u) - vpi(v);
            if vpi(x) > vpi(y)
                x = vpi(x) - vpi(y);
            else
                x = vpi(x) + vpi(p) - vpi(y);
            end
        else
            v = vpi(v) - vpi(u);
            if vpi(y) > vpi(x)
                y = vpi(y) - vpi(x);
            else
                y = vpi(y) + vpi(p) - vpi(x);
            end
        end
    end
    if u == 1
        R = mod(vpi(x), vpi(p));
%         R=x;
    elseif v == 1
        R = mod(vpi(y), vpi(p));
%         R=y;
    end
    
end