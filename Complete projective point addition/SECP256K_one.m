clear
clc
% get new points at https://andrea.corbellini.name/ecc/interactive/modk-add.html

%% Constants for SECP256K1

<<<<<<< HEAD
% primeP = 4999;
% XP = 1;
% YP = 2;
% b=3;
% m='2';

m = '2';
primeP = '115792089237316195423570985008687907853269984665640564039457584007908834671663';
XP = '55066263022277343669578718895168534326250603453777594175500187360389116729240';
YP = '32670510020758816978083085130507043184471273380659243275938904335757337482424';
b=7;
=======
primeP = 10001;
XP = 1;
YP = 2;
b=3;
m='2';

% m = '483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8';
% primeP = '115792089237316195423570985008687907853269984665640564039457584007908834671663';
% XP = '55066263022277343669578718895168534326250603453777594175500187360389116729240';
% YP = '32670510020758816978083085130507043184471273380659243275938904335757337482424';
% b=7;
>>>>>>> e6b93ed981e22d026e5264502b473722bbdeb0d8

ZP = 1;
%m_binary = hexToBinaryVector(m, [],'LSBFirst');
m_binary = flip(dec2bin(hex2dec(m)), 2);
i = length(m_binary)-1;
XR = XP;
YR = YP;
ZR = ZP;
b3 = b*3;
%% Using projective point addition exuations over a prime field 

tic; % start time
while i>0 
    [XR, YR, ZR] = projective_point_addition(XR,YR,ZR, XR,YR,ZR, b3, primeP);
    if m_binary(i)=='1'
        [XR, YR, ZR] = projective_point_addition(XR,YR,ZR, XP,YP,ZP, b3, primeP);
    end   
    i=i-1;
end

%% Perform modular multiplicative inverse
ZR
primeP
ZR_inverse = binary_inversion_algorithm(ZR, primeP)
x = mod(vpi(XR)*vpi(ZR_inverse), vpi(primeP))
y = mod(vpi(YR)*vpi(ZR_inverse), vpi(primeP))
simulation_time = toc; %stop time

fprintf('Simulation time: %f seocnds\n', simulation_time);



