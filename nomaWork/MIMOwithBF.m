%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All rights reserved by Krishna Pillai, http://www.dsplog.com
% The file may not be re-distributed without explicit authorization
% from Krishna Pillai.
% Checked for proper operation with Octave Version 3.0.0
% Author        : Krishna Pillai
% Email         : krishna@dsplog.com
% Version       : 1.0
% Date          : 12 April 2009
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script for computing the BER for BPSK modulation in a
% Rayleigh fading channel with and without transmit beamforming

clear
N = 10^6 % number of bits or symbols

% Transmitter
ip = rand(1,N)>0.5; % generating 0,1 with equal probability
s = 2*ip-1; % BPSK modulation 0 -> -1; 1 -> 0 
nTx = 2;
 
Eb_N0_dB = [0:35]; % multiple Eb/N0 values

for ii = 1:length(Eb_N0_dB)
   
   n = 1/sqrt(2)*[randn(1,N) + j*randn(1,N)]; % white gaussian noise, 0dB variance 
   h = 1/sqrt(2)*[randn(nTx,N) + j*randn(nTx,N)]; % Rayleigh channel

   sr = (1/sqrt(nTx))*kron(ones(nTx,1),s); 
   
   % Channel and noise Noise addition
   hEff = h.*exp(-j*angle(h));
   y1 = sum(h.*sr,1) + 10^(-Eb_N0_dB(ii)/20)*n; 
   y2 = sum(hEff.*sr,1) + 10^(-Eb_N0_dB(ii)/20)*n; 

   % equalization
   y1Hat = y1./sum(h,1); 
   y2Hat = y2./sum(hEff,1); 

   % receiver - hard decision decoding
   ip1Hat = real(y1Hat)>0;
   ip2Hat = real(y2Hat)>0;

   % counting the errors
   nErr1(ii) = size(find([ip- ip1Hat]),2);
   nErr2(ii) = size(find([ip- ip2Hat]),2);

end

simBer1 = nErr1/N; % simulated ber (no beam forming)
simBer2 = nErr2/N; % simulated ber (with beam forming)

theoryBerAWGN = 0.5*erfc(sqrt(10.^(Eb_N0_dB/10))); % theoretical ber
EbN0Lin = 10.^(Eb_N0_dB/10);
theoryBer = 0.5.*(1-sqrt(EbN0Lin./(EbN0Lin+1)));
p = 1/2 - 1/2*(1+1./EbN0Lin).^(-1/2);
theoryBer_nRx2 = p.^2.*(1+2*(1-p)); 


close all
figure
semilogy(Eb_N0_dB,theoryBer,'bp-','LineWidth',2);
hold on
semilogy(Eb_N0_dB,simBer1,'ys-','LineWidth',2);
semilogy(Eb_N0_dB,theoryBer_nRx2,'gp-','LineWidth',2);
semilogy(Eb_N0_dB,simBer2,'mx-','LineWidth',2);
axis([0 35 10^-5 0.5])
grid on
legend('1tx-1rx (theory)','2tx-1rx (no beamforming-sim)','1tx-2rx (mrc-theory)','2tx-1rx (beamforming-sim)');

xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('BER for BPSK modulation in Rayleigh channel');









