%code for OFDM signal transmission and reception in AWGN channel
% code developed by K.VINOTH BABU,VIT UNIVERSITY,VELLORE,INDIA on 07/10/2010
n = 256; % Number of bits to process
x = 2-randi(2,[n,1]); % Random binary data stream
x=x';
M = 16; % Size of signal constellation
k = log2(M); % Number of bits per symbol
xsym = bi2de(reshape(x,k,length(x)/k).','left-msb');% Convert the bits in x into k-bit symbols.
y = qammod(x,4,'gray'); % Modulate using 16-QAM.
tu=3.2e-6;%useful symbol period
tg=0.8e-6;%guard interval length
ts=tu+tg;%total symbol duration
nmin=0;
nmax=64;%total number of subcarriers
scb=312.5e3;%sub carrier spacing
fc=3.6e9;%carrier frequency
Rs=fc;
tt=0: 6.2500e-008:ts-6.2500e-008;
c=ifft(y,nmax);%IFFT 
s=real(c'.*(exp(1j*2*pi*fc*tt)));%bandpass modulation
figure;
plot(real(s),'b');title('OFDM signal transmitted');
figure;
plot(10*log10(abs(fft(s,nmax))));title('OFDM spectrum');
xlabel('frequency')
ylabel('power spectral density')
title('Transmit spectrum OFDM');
snr=10;%signal to noise ratio
ynoisy = awgn(s,snr,'measured');%awgn channel
figure;
plot(real(ynoisy),'b');title('received OFDM signal with noise');
z=ynoisy.*exp(1i*2*pi*fc*tt);%Bandpass demodulation
z=fft(z,nmax);%FFT
zsym = qamdemod(z,4,'gray');%demodulation of bandpass data
z = de2bi(zsym,'left-msb'); % Convert integers to bits.
z = reshape(z.',numel(z),1);%matrix to vector conversion
[noe,ber] = biterr(x,zsym) ;%BER calculation
figure;
subplot(211);stem(x(1:256));title('Original Message');
subplot(212);stem(z(1:256));title('recovered Message');
