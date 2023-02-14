clear
clc

primeP = '115792089237316195423570985008687907853269984665640564039457584007908834671663';
b = 7;
b3 = b*3;
XP = '55066263022277343669578718895168534326250603453777594175500187360389116729240';
YP = '32670510020758816978083085130507043184471273380659243275938904335757337482424';
ZP = 1;

m = '483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8';
m_binary = hexToBinaryVector(m, [],'LSBFirst');
i = length(m_binary)-1;
XR = XP;
YR = YP;
ZR = ZP;

tic;
while i>0
    
    [XR, YR, ZR] = projective_point_addition(XR,YR,ZR, XR,YR,ZR, b3, primeP);
    if m_binary(i)==1
        [XR, YR, ZR] = projective_point_addition(XR,YR,ZR, XP,YP,ZP, b3, primeP);
    end
    
    i=i-1;
end
simulation_time = toc;
%%
fprintf('Simulation time: %f seocnds\n', simulation_time);
XR
YR
ZR
%%





