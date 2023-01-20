function[y] = epsilon_1(x)
%epsilon 1
    a = rotr(x,6);
    b = rotr(x,11);
    d = rotr(x,25);
    
    c = x_or(a,b);
    
    y = x_or(c,d);   
end