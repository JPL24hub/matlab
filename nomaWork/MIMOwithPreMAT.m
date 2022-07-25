%MIMO with Zero Forcing equalizer
% Script for computing the BER for BPSK modulation in a
% Rayleigh fading channel with 2 Tx, 2Rx MIMO channel 
% Zero Forcing equalization

clear all
N_iter = 100; % number of bits or symbols
N=64;
Eb_N0_dB = 0:5:40; % multiple Eb/N0 values
nTx = 2;
nRx = 2;

for ii = 1:length(Eb_N0_dB)
    bittot=0;
    BER1=0;
    for j = 1:N_iter

    % Transmitter
    ip1 = rand(1,N)>0.5; % generating 0,1 with equal probability
    x1 = 2*ip1-1; % BPSK modulation 0 -> -1; 1 -> 0
    
    ip2 = rand(1,N)>0.5; % generating 0,1 with equal probability
    x2 = 2*ip2-1; % BPSK modulation 0 -> -1; 1 -> 0
    
    h11 = 1/sqrt(2)*(randn(1,N) + 1i*randn(1,N)); % Rayleigh channel
    h12 = 1/sqrt(2)*(randn(1,N) + 1i*randn(1,N)); % Rayleigh channel
    h21 = 1/sqrt(2)*(randn(1,N) + 1i*randn(1,N)); % Rayleigh channel
    h22 = 1/sqrt(2)*(randn(1,N) + 1i*randn(1,N)); % Rayleigh channel
    
    z1 = 1/sqrt(2)*(randn(1,N) + 1i*randn(1,N)); % white gaussian noise, 0dB variance
    z2 = 1/sqrt(2)*(randn(1,N) + 1i*randn(1,N)); % white gaussian noise, 0dB variance
    %[h11 h12
    % h21 h22]
    
       A = h11;
       B = h12;
       C = h21;
       D = h22;
       
       p1 = -D./(C.*B-D.*A);
       p3 = C./(C.*B-D.*A);

       p2 = -B./(D.*A-C.*B);
       p4 = A./(D.*A-C.*B);
       
       
       u1 = x1.*p1 + x2.*p2;
       u2 = x1.*p3 + x2.*p4;
       
    y1 = h11.*u1 + h12.*u2 + 10^(-Eb_N0_dB(ii)/20).*z1;
    y2 = h21.*u1 + h22.*u2 + 10^(-Eb_N0_dB(ii)/20).*z2;
       
    ipHat1 = real(y2)>0;
    BER1 = BER1 + size(find((ip1- ipHat1)),2);
    
    bittot = bittot + N;
    end

    % counting the errors
    simBer(ii) = BER1/bittot;
    
end

EbN0Lin = 10.^(Eb_N0_dB/10);
theoryBer_nRx1 = 0.5.*(1-1*(1+1./EbN0Lin).^(-0.5)); 
p = 1/2 - 1/2*(1+1./EbN0Lin).^(-1/2);
theoryBerMRC_nRx2 = p.^2.*(1+2*(1-p)); 
%%

figure
% semilogy(Eb_N0_dB,theoryBer_nRx1,'bp-','LineWidth',2)
% hold on
% semilogy(Eb_N0_dB,theoryBerMRC_nRx2,'kd-','LineWidth',2)
semilogy(Eb_N0_dB,simBer,'mo-','LineWidth',2)
% axis([0 25 10^-5 0.5])
grid on
legend('theory (nTx=1,nRx=1)', 'theory (nTx=1,nRx=2, MRC)', 'sim (nTx=2, nRx=2, ZF)')
xlabel('Average Eb/No,dB')
ylabel('Bit Error Rate')
title('BER for BPSK modulation with 2x2 MIMO and ZF equalizer (Rayleigh channel)')
% end




