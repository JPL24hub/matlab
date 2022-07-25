%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All rights reserved by Krishna Pillai, http://www.dsplog.com
% The file may not be re-distributed without explicit authorization
% from Krishna Pillai.
% Checked for proper operation with Octave Version 3.0.0
% Author        : Krishna Pillai
% Email         : krishna@dsplog.com
% Version       : 1.0
% Date          : 8 August 2008
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script for computing the BER for BPSK modulation in a
% Rayleigh fading channel

% clear all
N = 64; % number of bits or symbols
iter = 100;
Eb_N0_dB = [-3:35]; % multiple Eb/N0 values

for ii = 1:length(Eb_N0_dB)
    BER=0;
    bittot=0;
    for m=1:iter
    % Transmitter
        ip = rand(1,N)>0.5; % generating 0,1 with equal probability
        s = 2*ip-1; % BPSK modulation 0 -> -1; 1 -> 0 
   
        n = 1/sqrt(2)*[randn(1,N) + j*randn(1,N)]; % white gaussian noise, 0dB variance 
        h = 1/sqrt(2)*[randn(1,N) + j*randn(1,N)]; % Rayleigh channel
   
   % Channel and noise Noise addition
        y = diag(h)*diag(s) + 10^(-Eb_N0_dB(ii)/20)*diag(n); 

   % equalization
        yHat = y/diag(h);

   % receiver - hard decision decoding
        ipHat = real(diag(yHat))>0;

   % counting the errors
        BER = BER + size(find([ip- ipHat']),2);
        bittot = bittot + N;
    end
nErr(ii) = BER/bittot;
end

simBer = nErr/N; % simulated ber
theoryBerAWGN = 0.5*erfc(sqrt(10.^(Eb_N0_dB/10))); % theoretical ber
EbN0Lin = 10.^(Eb_N0_dB/10);
theoryBer = 0.5.*(1-sqrt(EbN0Lin./(EbN0Lin+1)));

%% plot
close all
figure
semilogy(Eb_N0_dB,theoryBerAWGN,'cd-','LineWidth',2);
hold on
semilogy(Eb_N0_dB,theoryBer,'bp-','LineWidth',2);
semilogy(Eb_N0_dB,nErr,'mx-','LineWidth',2);
axis([-3 35 10^-5 0.5])
grid on
legend('AWGN-Theory','Rayleigh-Theory', 'Rayleigh-Simulation');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('BER for BPSK modulation in Rayleigh channel');








