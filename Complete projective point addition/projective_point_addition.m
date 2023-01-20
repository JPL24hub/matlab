function[X3, Y3, Z3] = projective_point_addition(X1,Y1,Z1, X2,Y2,Z2, B3)
    t0 = X1*X2;
    t1 = Y1*Y2;
    t2 = Z1*Z2;
    t3 = X1 + Y1;
    t4 = X2 + Y2;
    t3 = t3*t4;
    
    
    
    X3 = 1;
    Y3 = 1;
    Z3 = 1;
    
end
