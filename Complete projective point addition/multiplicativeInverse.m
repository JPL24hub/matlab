function [inverse] = multiplicativeInverse(a, m)

    gcd_am = gcd(vpi(a), vpi(m)); %Find gcd

    if gcd_am ~= 1
        disp('Multiplicative invers does not exist');
        inverse = NaN;
    else
        [x, y] = extendedEuclidean(a, m) %fing x, y such that ax + my = gcd(a, m)=1.
        
        %The multiplicative inverse of a modulo m is given by x mod m
        %However, if x is negative, you need to add b to it to make it positive
        if (x<0)
            x = vpi(x)+vpi(m);
        end
        inverse = mod(vpi(x), vpi(m));
    end

end
