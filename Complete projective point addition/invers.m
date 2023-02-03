function[X, Y, Z] = invers(X, Y, Z, primeP)
Z_invers = Z^(primeP+2);
X = mod(vpi(X)*vpi(Z_invers), vpi(primeP));
Y = mod(vpi(Y)*vpi(Z_invers), vpi(primeP));
Z = mod(vpi(Z)*vpi(Z_invers), vpi(primeP));
end