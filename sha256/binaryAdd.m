function[sum, cary] = binaryAdd(a, b)
% adding binary
    i=length(a);
    
    [s,c] = bitAdd(a(1,i), b(1,i), 0);
    sum(1,i) = s;
    i=i-1;
    while (i>=1)
        [s, cary] = bitAdd(a(1,i), b(1,i), c);
        c = cary;
        sum(1,i) = s; 
        i=i-1;
    end
end