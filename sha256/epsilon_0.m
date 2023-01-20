function[y] = epsilon_0(x)
%epsilon 0
    a = rotr(x,2);
    b = rotr(x,13);
    d = rotr(x,22);
   
    c = x_or(a,b);
    
    y = x_or(c,d);   
end