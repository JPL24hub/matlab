function [g,a,b] = Extended_Euclidean_Int(v,u)
    r = Inf; 
    A = [v; u];
    q = [Inf]; % initalized with an used element. To match the number of rows
  
    while (r > 1) 
        c = floor(v/u);          
        r = v - c*u; 
        v = u;
        u = r;
        q = [q; c];  
        if r == 0
            break;
        end  
        A = [A; r];
    end
      
    [n,m] = size(A);
    y = zeros(n,1); 
    y(n-1: n,1) = [1; 0];    
    
    for i = n-1:-1:2
        y(i-1, 1) = q(i, 1) * y(i,1) + y(i+1, 1);
    end   
    
    r = A(n, 1);
    evaluatedSubtraction = A(1,1)*y(2,1) - A(2,1)*y(1,1);
    
    if r == evaluatedSubtraction
       b =  -1*y(1,1);
       a = y(2,1);
    else
       b = y(1,1);
       a = -1*y(2,1);
    end 
      
    g = r;
end