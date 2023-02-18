function[lcm] = test(m, n)
% least common multiple (LCM)
    a=m;
    b=n;
    while(a ~= b)
        if(a<b)
            a = a + m;
        else
            b = b + n;
        end
    end
    lcm = a;
end