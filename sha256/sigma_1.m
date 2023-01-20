function[y] = sigma_1(x)
%sigma 1
    a = rotr(x,17);
    b = rotr(x,19);
    d = shr(x,10);
   
    c = x_or(a,b);
    
    y = x_or(c,d);   
end