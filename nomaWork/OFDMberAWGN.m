% close all
% clear all
clc
nbitpersym  = 52;   % number of bits per OFDM symbol (same as the number of subcarriers for BPSK)
nsym        = 10^4; % number of symbols
len_fft     = 64;   % fft size
sub_car     = 52;   % number of data subcarriers
EbNo        = 0:5:40;
EsNo= EbNo + 10*log10(52/64)+ 10*log10(64/80); % symbol to noise ratio
snr= EsNo - 10*log10(64/80); % snr as to be used by awgn fn.

%%

% Generating data
t_data=2-randi(2,[nbitpersym*nsym,1]);

% modulating data
% mod_data = modulate(M,t_data);
mod_data = pskmod(t_data,2); % modulation object

% serial to parallel conversion
par_data = reshape(mod_data,nbitpersym,nsym).';
%%

% pilot insertion
pilot_ins_data=[zeros(nsym,6) par_data(:,[1:nbitpersym/2]) zeros(nsym,1) par_data(:,[nbitpersym/2+1:nbitpersym]) zeros(nsym,5)];

% fourier transform time doamain data and normalizing the data
IFFT_data = (64/sqrt(52))*ifft(fftshift(pilot_ins_data.')).';

% addition cyclic prefix
cylic_add_data = [IFFT_data(:,[49:64]) IFFT_data].';

% parallel to serial coversion
ser_data = reshape(cylic_add_data,80*nsym,1);
%%

% passing thru channel
h=rayleighchan(1/10000,10);
changain1=filter(h,ones(nsym*80,1));
a=max(max(abs(changain1)));
changain1=changain1./a;
chan_data = changain1.*ser_data;
no_of_error=[];
ratio=[];
%%
for ii=1:length(snr)
  
chan_awgn = awgn(chan_data,snr(ii),'measured'); % awgn addition

chan_awgn =a* chan_awgn./changain1;             % assuming ideal channel estimation

ser_to_para = reshape(chan_awgn,80,nsym).';     % serial to parallel coversion

cyclic_pre_rem = ser_to_para(:,[17:80]);        %cyclic prefix removal

FFT_recdata =(sqrt(52)/64)*fftshift(fft(cyclic_pre_rem.')).';                    % freq domain transform

%  FFT_recdata = FFT_recdata./FFT_recdata1;
rem_pilot = FFT_recdata (:,[6+[1:nbitpersym/2] 7+[nbitpersym/2+1:nbitpersym] ]); %pilot removal

ser_data_1 = reshape(rem_pilot.',nbitpersym*nsym,1);                             % serial coversion

%  ser_data_1  =  ser_data_1./abs(FFT_recdata1);
% z=modem.pskdemod(2); %demodulation object
demod_Data = pskdemod(ser_data_1,2);  %demodulating the data

[no_of_error(ii),ratio(ii)]=biterr(t_data,demod_Data) ; % error rate calculation
end
%%

% plotting the result
semilogy(EbNo,ratio,'--or','linewidth',2);
hold on;
EbN0Lin = 10.^(EbNo/10);
theoryBer = 0.5.*(1-sqrt(EbN0Lin./(EbN0Lin+1)));
semilogy(EbNo,theoryBer,'--ob','linewidth',2);
legend('simulated','theoritical')
grid on
xlabel('EbNo');
ylabel('BER')
title('Bit error probability curve for BPSK using OFDM');
