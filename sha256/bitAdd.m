function[sum, cary] = bitAdd(a, b, c)
% adding binary
    xr1 = xor(a,b);
    xr2 = xor(xr1,c);
    
    if a==1 && b==1 && c==1
        sum=1;
        cary=1;
    elseif a==0 && b==0 && c==0
        sum=0;
        cary=0;
    elseif xr2 == 0
        sum = 0;
        cary = 1;
    else
        sum = 1;
        cary = 0;
    end 
end