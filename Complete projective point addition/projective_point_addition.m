function[X3, Y3, Z3] = projective_point_addition(X1,Y1,Z1, X2,Y2,Z2, B3)
    t0 = X1 * X2;
    t1 = Y1 * Y2;
    t2 = Z1 * Z2;
    t3 = X1 + Y1;
    t4 = X2 + Y2;
    t3 = t3 * t4;
    t4 = t0 + t1;
    t3 = t3 - t4;
    t4 = Y1 + Z1;
    X3 = Y2 + Z2;
    t4 = t4 * X3;
    X3 = t1 + t2;
    t4 = t4 * X3;
    X3 = X1 + Z1;
    Y3 = X2 + Z2;
    X3 = X3 * Y3;
    Y3 = t0 + t2;
    Y3 = X3 * Y3;
    X3 = t0 + t0;
    t0 = X3 + t0;
    t2 = B3 * t2;
    Z3 = t1 + t2;
    t1 = t1 * t2;
    Y3 = B3 * Y3;
    X3 = t4 * Y3;
    t2 = t3 * t1;
    X3 = t2 * X3;
    Y3 = Y3 * t0;
    t1 = t1 * Z3;
    Y3 = t1 + Y3;
    t0 = t0 * t3;
    Z3 = Z3 * t4;
    Z3 = Z3 + t0;
    
    
end
