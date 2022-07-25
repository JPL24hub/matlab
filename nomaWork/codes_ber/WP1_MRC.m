
function [F1wp1,X1wp1,F2wp1,X2wp1, BERu11] = WP1_MRC(iter)

warning off;
clc

Nbps=1; 
M=2^Nbps;           % Modulation order=2/4/6 for QPSK/16QAM/64QAM
EbN0=0:5:40;        % EbN0, SNR
N_iter=iter;       % no of packets   % Number of iterations for each EbN0
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
     BERu1=0;
     bittot=0;
     
     for m=1:N_iter
         %% Round 1 antenna 1
          %h_ab  a=usr b=antena
      %user1
       channel=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);
       h11=zeros(1,Lch); 
       h11(Delay+1)=channel; % cir: channel impulse response
       h11r1=h11;
       H11r1=fft([h11r1 zeros(1,Nfft-Lch)]); % Channel frequency response
      %user2
       channel=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);
       h21=zeros(1,Lch); 
       h21(Delay+1)=channel; % cir: channel impulse response
       h21r1=h21;
       H21r1=fft([h21r1 zeros(1,Nfft-Lch)]); % Channel frequency response
      
       %Eve
       channel=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);
       hE1=zeros(1,Lch); 
       hE1(Delay+1)=channel; % cir: channel impulse response
       hEr1=hE1;
       HEr1=fft([hEr1 zeros(1,Nfft-Lch)]); % Channel frequency response
       %% Round 2 antenna 2
       %h_ab  a=usr b=antena
   
      %user1  
       channel=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);	
       h12=zeros(1,Lch); 
       h12(Delay+1)=channel; % cir: channel impulse response
       h12r2=h12;
       H12r2=fft([h12r2 zeros(1,Nfft-Lch)]); % Channel frequency response
      
      % user2
       channel=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);	
       h22=zeros(1,Lch); 
       h22(Delay+1)=channel; % cir: channel impulse response
       h22r2=h22;
       H22r2=fft([h22r2 zeros(1,Nfft-Lch)]); % Channel frequency response
     
       %Eve
       channel=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);	
       hE2=zeros(1,Lch); 
       hE2(Delay+1)=channel; % cir: channel impulse response
       hEr2=hE2;
       HEr2=fft([hEr2 zeros(1,Nfft-Lch)]); % Channel frequency response
       
       %% Generate Signal for each user/tx, modulate, and add CP
       
      % User1 Tx______________________________________________________________
      X= randi([0,M-1], Nused*Nframe,1)'; % bit: integer vector (1)
      Xmod= qammod(X,M,'gray')/norms(Nbps); 	
      
      % User2 Tx______________________________________________________________
      X2= randi([0,M-1], Nused*Nframe,1)'; % bit: integer vector
      Xmod2= qammod(X2,M,'gray')/norms(Nbps);
      

      %%Modulted signal for user1 and user2      
      U1=Xmod;
      U2=Xmod2;
      %% precoder matrices
      %h_ab  a=usr b=antena

      M_1a=abs(diag(H22r2)).^(2)*(abs(diag(H22r2)).^(2)*abs(diag(H11r1)).^(2)-abs(diag(H21r1)).^(2)*abs(diag(H12r2)).^(2))^-1;
      M_1b=-(M_1a*abs(diag(H21r1)).^(2))/abs(diag(H22r2).^(2));

      M_2a=abs(diag(H12r2)).^(2)*(abs(diag(H21r1)).^(2)*abs(diag(H12r2)).^(2)-abs(diag(H22r2)).^(2)*abs(diag(H11r1)).^(2))^-1;
      M_2b=-(M_2a*abs(diag(H11r1)).^(2))/abs(diag(H12r2)).^(2);
           
       %% Super position
        %%Antenna 1
      Tx1=sqrt(1)*U1*M_1a+sqrt(1)*U2*M_2a;
        %Antenna 2    
      Tx2=sqrt(1)*U1*M_1b+sqrt(1)*U2*M_2b;
      
     %% FOR PAPR
     
     x1np= ifft(U1);
     x_GI1np= guard_interval(Ng,Nfft,NgType,x1np);       
     x2np= ifft(U2);
     x_GI2np= guard_interval(Ng,Nfft,NgType,x2np); 

       %% 
       %Tx1
         x1= ifft(Tx1)/norm(ifft(Tx1),2);  %% This is how i am normalizign it :step1
         norm(ifft(Tx1),2);
        % pause()
         x_GI1= guard_interval(Ng,Nfft,NgType,x1);       
         y11= conv(x_GI1,h11r1); %u1 r1
         
         %EbN0=20;
         snr = EbN0(k)+10*log10(Nbps*(Nused/Nfft)); % SNR vs. Eb/N0
        %snr=snr/2;
         EEs=sum(x1).^2/length(x1);    
         %EEs=1;
         noise_mag  = sqrt(0.5*EEs*((10.^(-snr/10))));
         y_G11r1 = norm(ifft(Tx1),2)*y11+ noise_mag*(randn(size(y11))+1i*randn(size(y11))); %% This is how i am normalizign : step2
         Y11= fft(remove_GI(Ng,Nsym,NgType,y_G11r1));
         %%
       %Tx2 
         x2= ifft(Tx2)/norm(ifft(Tx2),2);
         x_GI2= guard_interval(Ng,Nfft,NgType,x2);       
         y12=conv(x_GI2,h12r2); %u1 r2
  
        EEs=sum(x2).^2/length(x2);    
        noise_mag  = sqrt(0.5*EEs*((10.^(-snr/10))));  
        y_G12r2 = norm(ifft(Tx2),2)* y12+ noise_mag*(randn(size(y12))+1i*randn(size(y12)));
 
        Y12= fft(remove_GI(Ng,Nsym,NgType,y_G12r2));
         %%
         
         
                  % %%%%%%% FOR PAPR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         if k==1
         PAPRtx1(1,m)=max((abs(x_GI1)).^2)/mean((abs(x_GI1)).^2);
         PAPRtx2(1,m)=max((abs(x_GI2)).^2)/mean((abs(x_GI2)).^2);
         PAPRtx3(1,m)=max((abs(x_GI1np)).^2)/mean((abs(x_GI1np)).^2);
         PAPRtx4(1,m)=max((abs(x_GI2np)).^2)/mean((abs(x_GI2np)).^2);
         %pause();
         end
         
         % %%%%%%% FOR PAPR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
         
         
         YU1= (H11r1)'.*transpose(Y11)+(H12r2)'.*transpose(Y12);
         
         %%Demodulate
         UU1i=qamdemod(YU1*norms(Nbps),M,'gray');

          %% BER
          
          BERu1= BERu1+sum(sum(de2bi(UU1i,Nbps)~=de2bi(X,Nbps)));         
          bittot=bittot+Nfft*Nbps; 
         
     end
     
     %BER
      BERu11(k)=BERu1/bittot;
     
      
end
%%
% %BER
% figure(1);
% semilogy(EbN0,BERu11,'or--','LineWidth',2);
% xlabel('SNR (dB)');
% ylabel('Bit Error Rate(BER)');
% grid on;
% legend('BER(user-1)');


%PAPR
PAPRtx1dB=10*log10(PAPRtx1); % in dB
[F1wp1,X1wp1] = ecdf(PAPRtx1dB); %CDF

PAPRtx2dB=10*log10(PAPRtx2); % in dB
[F2wp1,X2wp1] = ecdf(PAPRtx2dB); %CDF



end