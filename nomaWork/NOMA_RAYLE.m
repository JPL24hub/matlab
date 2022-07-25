
%  function [ber_theo1,ber_theo2] = NOMA_BER_Rayleigh()

%ber_th1 is far user, ber_th2 is near user
%https://ecewireless.blogspot.com/2020/04/how-to-simulate-ber-capacity-and-outage.html
clc; 
clear all

N = 10^3;

d1 = 1000; d2 = 500;    %Distances of users from base station (BS)
a1 = 0.75; a2 = 0.25;   %Power allocation factors
eta = 4;                %Path loss exponent

%Generate rayleigh fading coefficient for both users
h11 = sqrt(d1^-eta)*(randn(1,N)+1i*randn(1,N))/sqrt(2);
h21 = sqrt(d2^-eta)*(randn(1,N)+1i*randn(1,N))/sqrt(2);

h12 = sqrt(d1^-eta)*(randn(1,N)+1i*randn(1,N))/sqrt(2);
h22 = sqrt(d2^-eta)*(randn(1,N)+1i*randn(1,N))/sqrt(2);


% EbN0=0:5:40;        % EbN0, SNR
Pt = 0:5:40;                %Transmit power in dBm
pt = (10^-3)*10.^(Pt/10);   %Transmit power in linear scale
BW = 10^6;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

%Generate noise samples for both users
w11 = sqrt(no)*(randn(1,N)+1i*randn(1,N))/sqrt(2);
w21 = sqrt(no)*(randn(1,N)+1i*randn(1,N))/sqrt(2);

w12 = sqrt(no)*(randn(1,N)+1i*randn(1,N))/sqrt(2);
w22 = sqrt(no)*(randn(1,N)+1i*randn(1,N))/sqrt(2);

%Generate random binary data for two users
data1 = randi([0 1],1,N);  %Data bits of user 1
data2 = randi([0 1],1,N);  %Data bits of user 2

%Do BPSK modulation of data
x1 = 2*data1 - 1;
x2 = 2*data2 - 1;

p = length(Pt);
for u = 1:p
    %Do superposition coding
    
    A = diag(h11);
    B = diag(h12);
    C = diag(h21);
    D = diag(h22);
        
    S = ((A+B)*diag(x2));
    T = ((C+D)*diag(x1));
        
        
    r2 = ((S*C)-(A*T))/((D*A)-(C*B));
    r1 = (-S-(B*r2))/A;
    
    Tx1 = sqrt(pt(u))*(sqrt(a1)*x1 + sqrt(a2)*x2 + r1);
    Tx2 = sqrt(pt(u))*(sqrt(a1)*x1 + sqrt(a2)*x2 + r2);
    
    %Received signals
    y11 = h11.*Tx1 + w11;
    y21 = h21.*Tx1 + w21;
    y12 = h12.*Tx2 + w12;
    y22 = h22.*Tx2 + w22;
    
    y1 = y11 + y12;
    y2 = y21 + y22;
    
    %Equalize 
    eq1 = y11./(h11 + h12);
    eq1 = diag(eq1);
    eq1 = eq1';
    eq2 = y21./(h21 + h22);
    eq2 = diag(eq2);
    eq2 = eq2';
   
    
    %AT USER 1--------------------
    %Direct decoding of x1 from y1
    x1_hat = zeros(1,N);
    x1_hat(eq1>0) = 1;

    
    %Compare decoded x1_hat with data1 to estimate BER
   
    ber1(u) = biterr(data1,x1_hat)/N; % biterr returns howmany bits are not correctly received
 
    %----------------------------------
    
    %AT USER 2-------------------------
    %Direct decoding of x1 from y2
    x12_hat = ones(1,N);
    x12_hat(eq2<0) = -1;
    
    
    y2_dash = eq2 - sqrt(a1*pt(u))*x12_hat;  
    x2_hat = zeros(1,N);
    x2_hat(real(y2_dash)>0) = 1;

    ber2(u) = biterr(x2_hat, data2)/N;
    %-----------------------------------   
    
  
end
%%
% SNR = Pt/No
%%
semilogy(Pt, ber1,'r', 'linewidth',1.5); hold on; grid on;
semilogy(Pt, ber2,'b', 'linewidth',1.5);
xlabel('Transmit power (P in dBm)');
ylabel('BER');
legend('Sim. User 1/Far user','Sim. User 2/Near user','Theo. User 1/Far user','Theo. User 2/Near user');
%   end











