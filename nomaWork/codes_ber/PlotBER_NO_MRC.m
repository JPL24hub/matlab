
clear all;
warning off;
%%

Nbps=1; 
M=2^Nbps;  % Modulation order=2/4/6 for QPSK/16QAM/64QAM
EbN0=[0:5:40];    % EbN0, SNR
%CR=1;%1/2;%CODING RATE
%EbN0=[40];   % EbN0
%se=1; %with MRc/ without MRc
%
N_iter=4000%1e5;   % no of packets   % Number of iterations for each EbN0
Nframe=1;         % no. of sybmols in on % Number of symbols per frame
sigPow=0;         % Signal power initialization

norms=[1 sqrt(2) 0 sqrt(10) 0 sqrt(42)]; 

%% Bob
PowerdB=[0 -8 -17 -21 -25]; % Channel tap power profile 'dB'
Delay=[0 3 5 6 8];          % Channel delay 'sample'
Power=10.^(PowerdB/10);     % Channel tap power profile 'linear scale'
Power=Power/(sum(Power));
Ntap=length(PowerdB);       % Chanel tap number
Lch=Delay(end)+1;           %Channel length

%% Eve
PowerdBE=[0 -8 -17 -21 -25];
DelayE=[0 3 5 6 8];
PowerE=10.^(PowerdBE/10);     % Channel tap power profile 'linear scale'
NtapE=length(PowerdBE);  % Chanel tap number
LchE=DelayE(end)+1;      % Channel length(may be diff than no taps)
%%

Nfft=64;           % FFT size
Ng=16; %Ng=3; %    % Ng=0: Guard interval length
NgE=16;
%Ng=0;%Nfft/4;
Nsym=Nfft+Ng;      % Symbol duration
Nvc=0;%Nfft/4;        % Nvc=0: no virtual carrier
Nused=Nfft ; %-Nvc;

NgType=1; % NgType=1/2 for cyclic prefix/zero padding
if NgType==1, nt='CP';  elseif NgType==2, nt='ZP';   end
Ch=1;  % Ch=0/1 for AWGN/multipath channel
if Ch==0, chType='AWGN'; Target_neb=100; else chType='CH'; Target_neb=500; end
%figure(Ch+1), clf

%%

for k=1:length(EbN0)
     BERu1=0;
     BERu2=0;
     
     BERI1E=0;
     BERI2E=0;
     
     BERu1E=0;
     BERu2E=0;
     
     bittot=0;
     
     for m=1:N_iter
         %% Round 1 antenna 1
          %h_ab  a=usr b=antena
      %user1
       channel=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);%%%%%
       h11=zeros(1,Lch); 
       h11(Delay+1)=channel; % cir: channel impulse response
       h11r1=h11;
       H11r1=fft([h11r1 zeros(1,Nfft-Lch)]); % Channel frequency response
      %user2
       channel=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);%%%%%
       h21=zeros(1,Lch); 
       h21(Delay+1)=channel; % cir: channel impulse response
       h21r1=h21;
       H21r1=fft([h21r1 zeros(1,Nfft-Lch)]); % Channel frequency response
      
       %Eve
       channel=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);%%%%%
       hE1=zeros(1,Lch); 
       hE1(Delay+1)=channel; % cir: channel impulse response
       hEr1=hE1;
       HEr1=fft([hEr1 zeros(1,Nfft-Lch)]); % Channel frequency response
       %% Round 2 antenna 2
       %h_ab  a=usr b=antena
   
      %user1  
       channel=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);%%%%%
       h12=zeros(1,Lch); 
       h12(Delay+1)=channel; % cir: channel impulse response
       h12r2=h12;
       H12r2=fft([h12r2 zeros(1,Nfft-Lch)]); % Channel frequency response
      
      % user2
       channel=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);%%%%%
       h22=zeros(1,Lch); 
       h22(Delay+1)=channel; % cir: channel impulse response
       h22r2=h22;
       H22r2=fft([h22r2 zeros(1,Nfft-Lch)]); % Channel frequency response
     
       %Eve
       channel=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);%%%%%
       hE2=zeros(1,Lch); 
       hE2(Delay+1)=channel; % cir: channel impulse response
       hEr2=hE2;
       HEr2=fft([hEr2 zeros(1,Nfft-Lch)]); % Channel frequency response
       %% Generate Signal for each user/tx, modulate, and add CP
       
      % User1 Tx______________________________________________________________
      X= randi([0,M-1], Nused*Nframe,1)'; % bit: integer vector (1)
      %   X= randint(1,Nused*Nframe,M); % bit: integer vector %Non binary data 0 1 2 3
      Xmod= qammod(X,M,0,'gray')/norms(Nbps); %% (2)
      if NgType~=2, x_GI=zeros(1,Nframe*Nsym); %(3)
      elseif NgType==2, x_GI= zeros(1,Nframe*Nsym+Ng);
        % Extend an OFDM symbol by Ng zeros 
      end
      
      % User2 Tx______________________________________________________________
      X2= randi([0,M-1], Nused*Nframe,1)'; % bit: integer vector
      %   X= randint(1,Nused*Nframe,M); % bit: integer vector %Non binary data 0 1 2 3
      Xmod2= qammod(X2,M,0,'gray')/norms(Nbps); %%
      if NgType~=2, x_GI2=zeros(1,Nframe*Nsym);
      elseif NgType==2, x_GI2= zeros(1,Nframe*Nsym+Ng);
        % Extend an OFDM symbol by Ng zeros 
      end
      %%Modulted signal for user1 and user2      
      U1=Xmod;
      U2=Xmod2;
      %% Precoding matrix 
      %h_ab  a=usr b=antena
      
       M_1a=diag(H22r2)*(diag(H22r2)*diag(H11r1)-diag(H21r1)*diag(H12r2))^-1;
       M_1b=-(M_1a*diag(H21r1))/diag(H22r2);
       
       M_2a=diag(H12r2)*(diag(H21r1)*diag(H12r2)-diag(H22r2)*diag(H11r1))^-1;
       M_2b=-(M_2a*diag(H11r1))/diag(H12r2);
       
        %% Normaization factore
       snr = EbN0(k)+10*log10(Nbps*(Nused/Nfft)); % SNR vs. Eb/N0
       noise_var = Nfft*0.5*10^(-snr/10); 
       
        %% Super position
        %%Antenna 1
       Tx1=sqrt(0.5)*U1*M_1a+sqrt(0.5)*U2*M_2a;
        %%Antenna 2    
       Tx2=sqrt(0.5)*U1*M_1b+sqrt(0.5)*U2*M_2b;
       
       
       %% PRecoded: This is imp (from X to Y)
       
       %Tx1
         x1= ifft(Tx1)/norm(ifft(Tx1),2);  %% This is how i am normalizign it :step1
         norm(ifft(Tx1),2);
        % pause()
         x_GI1= guard_interval(Ng,Nfft,NgType,x1);       
         y11= conv(x_GI1,h11r1); %u1 r1
         y21=(conv(x_GI1,h21r1)); %u2 r1
         yE1= conv(x_GI1,hEr1); %u2 r1
         %EbN0=20;
         snr = EbN0(k)+10*log10(Nbps*(Nused/Nfft)); % SNR vs. Eb/N0
        %snr=snr/2;
         EEs=sum(x1).^2/length(x1);    
         %EEs=1;
         noise_mag  = sqrt(0.5*EEs*((10.^(-snr/10))));
         y_G11r1 = norm(ifft(Tx1),2)*y11+ noise_mag*(randn(size(y11))+1i*randn(size(y11))); %% This is how i am normalizign : step2
         y_G21r1 = norm(ifft(Tx1),2)*y21+ noise_mag*(randn(size(y21))+1i*randn(size(y21)));
         y_G21E1 = norm(ifft(Tx1),2)*yE1+ noise_mag*(randn(size(yE1))+1i*randn(size(yE1)));
 
         
         Y11= fft(remove_GI(Ng,Nsym,NgType,y_G11r1));
         Y21= fft(remove_GI(Ng,Nsym,NgType,y_G21r1));
         YE1= fft(remove_GI(Ng,Nsym,NgType,y_G21E1));
         %%
       %Tx2 
         x2= ifft(Tx2)/norm(ifft(Tx2),2);
         x_GI2= guard_interval(Ng,Nfft,NgType,x2);       
         y12=conv(x_GI2,h12r2); %u1 r2
         y22= conv(x_GI2,h22r2); %u2 r2
         yE2= conv(x_GI2,hEr2); %u2 r1    
 
         %EbN0=10;
  
        EEs=sum(x2).^2/length(x2);    
        noise_mag  = sqrt(0.5*EEs*((10.^(-snr/10))));  
        y_G12r2 = norm(ifft(Tx2),2)* y12+ noise_mag*(randn(size(y12))+1i*randn(size(y12)));
        y_G22r2 = norm(ifft(Tx2),2)* y22+ noise_mag*(randn(size(y22))+1i*randn(size(y22)));
        y_G21E2 = norm(ifft(Tx2),2)*yE2+ noise_mag*(randn(size(yE2))+1i*randn(size(yE2)));
 
        Y12= fft(remove_GI(Ng,Nsym,NgType,y_G12r2));
        Y22= fft(remove_GI(Ng,Nsym,NgType,y_G22r2));
        YE2= fft(remove_GI(Ng,Nsym,NgType,y_G21E2));
         
         %%
         YU1= (Y11)+(Y12);
         YU2= (Y21)+(Y22);
         YUe= YE1+YE2;
         
         %%Demodulate
         UU1i=qamdemod(YU1*norms(Nbps),M,0,'gray');
         UU2i=qamdemod(YU2*norms(Nbps),M,0,'gray');
         UUee=qamdemod(YUe*norms(Nbps),M,0,'gray');
         
         %% Compare Output
          user1n=length(find(UU1i==X));
          user2n=length(find(UU2i==X2));
          
            %%External Eve
          user1en=length(find(UUee==X));
          user2en=length(find(UUee==X2));
            %%Internal Eve
          userI1en=length(find(X2==X));
          userI2en=length(find(X==X2));
          
          %% BER
          
          BERu1=BERu1+sum(sum(de2bi(UU1i,Nbps)~=de2bi(X,Nbps)));
          BERu2= BERu2+sum(sum(de2bi(UU2i,Nbps)~=de2bi(X2,Nbps)));
          
          
          %%External Eve
          BERu1E=BERu1E+sum(sum(de2bi(UUee,Nbps)~=de2bi(X,Nbps)));
          BERu2E=BERu2E+sum(sum(de2bi(UUee,Nbps)~=de2bi(X2,Nbps)));
 
          %%Internal Eve
          BERI1E=BERI1E+sum(sum(de2bi(X,Nbps)~=de2bi(X2,Nbps)));
          BERI2E=BERI2E+sum(sum(de2bi(X2,Nbps)~=de2bi(X,Nbps)));
          
          bittot=bittot+Nfft*Nbps; 
         
     end
     
     %BER
      BERu11(k)=BERu1/bittot;
      BERu21(k)=BERu2/bittot;
      BERu1E1(k)=BERu1E/bittot;
      BERu2E1(k)=BERu2E/bittot;
      
end

figure(1);
semilogy(EbN0,BERu11,'or-','LineWidth',2); hold on;
semilogy(EbN0,BERu21,'*b-','LineWidth',2);hold on;
semilogy(EbN0,BERu1E1,'vc-','LineWidth',2); hold on;
semilogy(EbN0,BERu2E1,'^r-','LineWidth',2); 
xlabel('SNR (dB)');
ylabel('Bit Error Rate(BER)');
legend('BER(Bob1)','BER(Bob2)','BER (E-Eve-1)','BER (E-Eve-2)');
