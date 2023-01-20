function[y] = sigma_0(x)
%sigma zero
    a = rotr(x,7);
    b = rotr(x,18);
    d = shr(x,3);
    
    c = x_or(a,b);
    y = x_or(c,d);   
end