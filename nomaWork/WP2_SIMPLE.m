
warning off;
clc

%%
Nbps=1; 
M=2^Nbps;           % Modulation order=2/4/6 for QPSK/16QAM/64QAM
EbN0=0:5:25;        % EbN0, SNR
N_iter=1000;       % no of packets   % Number of iterations for each EbN0
Nframe=1;           % no. of sybmols in on % Number of symbols per frame

norms=[1 sqrt(2) 0 sqrt(10) 0 sqrt(42)]; 

%% Bob,Eve
PowerdB=[0 -8 -17 -21 -25]; % Channel tap power profile 'dB'
Delay=[0 3 5 6 8];          % Channel delay 'sample'
Power=10.^(PowerdB/10);     % Channel tap power profile 'linear scale'
Power=Power/(sum(Power));
Ntap=length(PowerdB);       % Chanel tap number
Lch=Delay(end)+1;           %Channel length
%%

Nfft=64;           % FFT size
Ng=16;             % Ng=0: Guard interval length
Nsym=Nfft+Ng;      % Symbol duration
Nused=Nfft ;

NgType=1; % NgType=1/2 for cyclic prefix/zero padding
if NgType==1, nt='CP';  elseif NgType==2, nt='ZP';   end


%%

for k=1:length(EbN0)

     BER_B=0;
     BER_E=0;
     bittot=0;    
     
     for m=1:N_iter
         %% CHANNEL
      %HB1
       channel=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);
       hb1=zeros(1,Lch); 
       hb1(Delay+1)=channel; % cir: channel impulse response
       hb1r1=hb1;
       HB1=fft([hb1r1 zeros(1,Nfft-Lch)]); % Channel frequency response
       
      %HB2
       channel=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);
       nb2=zeros(1,Lch); 
       nb2(Delay+1)=channel; % cir: channel impulse response
       hb2r1=nb2;
       HB2=fft([hb2r1 zeros(1,Nfft-Lch)]); % Channel frequency response
    
       %HE1
       channel=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);
       he1=zeros(1,Lch); 
       he1(Delay+1)=channel; % cir: channel impulse response
       he1r1=he1;
       HE1=fft([he1r1 zeros(1,Nfft-Lch)]); % Channel frequency response
 
      %HE2 
       channel=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);	
       he2=zeros(1,Lch); 
       he2(Delay+1)=channel; % cir: channel impulse response
       he2r2=he2;
       HE2=fft([he2r2 zeros(1,Nfft-Lch)]); % Channel frequency response
       
       %%
       
      % User1 Tx______________________________________________________________
      xrand = randi([0,M-1], Nused*Nframe,1)'; % bit: integer vector (1)
      X = qammod(xrand,M,'gray')/norms(Nbps); 	
    

      %% Auxiliary signals
      ANP = 10^(-(EbN0(k))/10);
      W = randi([0 1], 1, 64);
      Q = randi([0 1], 1, 64);
      
      An = sqrt(ANP/2)*((2*W-1) + 1i*(2*Q-1));
      
      R1 = An/HB1;
      R2 = An/HB2;
      
      X1 = X + R1;
      X2 = X - R2;

       %% 
       %Tx1
         x1= ifft(X1)/norm(ifft(X1),2);  %% This is how i am normalizign it :step1
         norm(ifft(X1),2);
         
        % pause()
         x_GI1= guard_interval(Ng,Nfft,NgType,x1);       
         y1b= conv(x_GI1,hb1r1); %u1 r1
         y1e= conv(x_GI1,he1r1); %u2 r1
         
         %EbN0=20;
         snr = EbN0(k)+10*log10(Nbps*(Nused/Nfft)); % SNR vs. Eb/N0
        %snr=snr/2;
         EEs=sum(x1).^2/length(x1);    
         %EEs=1;
         noise_mag  = sqrt(0.5*EEs*((10.^(-snr/10))));
         
         
         y_G11r1 = norm(ifft(X1),2)*y1b + noise_mag*(randn(size(y1b))+1i*randn(size(y1b))); %% This is how i am normalizign : step2
         y_G21E1 = norm(ifft(X1),2)*y1e + noise_mag*(randn(size(y1e))+1i*randn(size(y1e)));
         
         Y1B= fft(remove_GI(Ng,Nsym,NgType,y_G11r1));
         Y1E= fft(remove_GI(Ng,Nsym,NgType,y_G21E1));
         %%
       %Tx2 
         x2= ifft(X2)/norm(ifft(X2),2);
         x_GI2= guard_interval(Ng,Nfft,NgType,x2);
         
         y2b = conv(x_GI2,hb2r1); %u1 r2
         y2e = conv(x_GI2,he2r2); %u2 r1 
 
         %EbN0=10;
        EEs=sum(x2).^2/length(x2);    
        noise_mag  = sqrt(0.5*EEs*((10.^(-snr/10)))); 

        y_G12r2 = norm(ifft(X2),2)*y2b + noise_mag*(randn(size(y2b))+1i*randn(size(y2b)));
        y_G21E2 = norm(ifft(X2),2)*y2e + noise_mag*(randn(size(y2e))+1i*randn(size(y2e)));
 
        Y2B= fft(remove_GI(Ng,Nsym,NgType,y_G12r2));
        Y2E= fft(remove_GI(Ng,Nsym,NgType,y_G21E2));
         
         %%
         YB= Y1B + Y2B;
         YE= Y1E + Y2E; 
         
         %%Demodulate
         UU1i=qamdemod(YB*norms(Nbps),M,'gray');
         UUee=qamdemod(YE*norms(Nbps),M,'gray');    
          %% BER
          
          BER_B = BER_B+sum(sum(de2bi(UU1i,Nbps)~=de2bi(xrand,Nbps)));    
          %%Internal Eve
          BER_E = BER_E+sum(sum(de2bi(UUee,Nbps)~=de2bi(xrand,Nbps))); %data send to user-2 seen by user-1   
          
          bittot=bittot+Nfft*Nbps; 
         
     end
     
     
     %BER
      BERu11(k) = BER_B/bittot;
      BERIu11(k) = BER_E/bittot;
      
end
%%

%BER
figure(1);
semilogy(EbN0,BERu11,'or--','LineWidth',2);  hold on;
semilogy(EbN0,BERIu11,'g--','LineWidth',2); hold on; %user-1 is eavs
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
grid on;
legend('user-1','user-2');





