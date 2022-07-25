

% function [BERu11, BERu21] = WP3_NO_MRC(iter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Version1: OFDM-With gprecoding
%%Updated: No channel at RX: NO MRC

%% Normalization method
 % x1= ifft(Tx1)/norm(ifft(Tx1),2); 
 %y_G11r1 = norm(ifft(Tx1),2)*y11+ noise_mag*(randn(size(y11))+j*randn(size(y11)));
  
warning off;
clc
clear all

%% ARQ parameters0
%Important parameters: L , Modulation, chanel coding , fft size, frmae size 
%L=2; %No of retransmsiion

Nbps=1; 
M=2^Nbps;  % Modulation order=2/4/6 for QPSK/16QAM/64QAM
EbN0=1:1:30;   % EbN0
%CR=1;%1/2;%CODING RATE
%EbN0=[40];   % EbN0
%se=1; %with MRc/ without MRc
%
N_iter=4000;%1e5;   % no of packets   % Number of iterations for each EbN0
Nframe=1;         % no. of sybmols in on % Number of symbols per frame
sigPow=0;         % Signal power initialization
%file_name=['OFDM_BER_' chType '_' nt '_' 'GL' num2str(Ng) '.dat'];
%fid=fopen(file_name, 'w+');
norms=[1 sqrt(2) 0 sqrt(10) 0 sqrt(42)];     % BPSK 4-QAM 16-QAM
%
NgType=1; % NgType=1/2 for cyclic prefix/zero padding
if NgType==1, nt='CP';  elseif NgType==2, nt='ZP';   end
Ch=1;  % Ch=0/1 for AWGN/multipath channel
if Ch==0, chType='AWGN'; Target_neb=100; else chType='CH'; Target_neb=500; end
%figure(Ch+1), clf

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

Nfft=64;           % FFT size
Ng=16; %Ng=3; %    % Ng=0: Guard interval length
NgE=16;
%Ng=0;%Nfft/4;
Nsym=Nfft+Ng;      % Symbol duration
Nvc=0;%Nfft/4;        % Nvc=0: no virtual carrier
Nused=Nfft ; %-Nvc;
 
for i=1:length(EbN0)
    
SERu1=0;
SERu2=0;
SERu1E=0;
SERu2E=0;

SERuI1E=0;
SERuI2E=0;

symtot=0;
    %Neb=0; Ntb=0; % Initialize the number of error/total bits
 BERu1=0;
 BERu2=0;
 
 BERI1E=0;
 BERI2E=0;

 BERu1E=0;
 BERu2E=0;
 
 
     TERu1=0;
     TERu2=0;
     TERu1E=0;
     TERu2E=0;
     TERI1E=0;
     TERI2E=0;
 
 bittot=0;
 
   for m=1:N_iter
       
       %%%
       s=0;% reset to begining  
       Packet_err_i_b=1;

       %h_ab  a=usr b=antena  r1 round 1
     %% Round 1
      %user1-atenna1 
        kk4=Nused/2+Nvc+1:Nfft; kk5=(Nvc~=0)+[1:Nused/2];
 
       channel=(randn(1,Ntap)+j*randn(1,Ntap)).*sqrt(Power/2);%%%%%
       h11=zeros(1,Lch); 
       h11(Delay+1)=channel; % cir: channel impulse response
       h11r1=h11;
       H11r1=fft([h11r1 zeros(1,Nfft-Lch)]); % Channel frequency response
      %user2
       channel=(randn(1,Ntap)+j*randn(1,Ntap)).*sqrt(Power/2);%%%%%
       h21=zeros(1,Lch); 
       h21(Delay+1)=channel; % cir: channel impulse response
       h21r1=h21;
       H21r1=fft([h21r1 zeros(1,Nfft-Lch)]); % Channel frequency response
      
       %Eve
       channel=(randn(1,Ntap)+j*randn(1,Ntap)).*sqrt(Power/2);%%%%%
       hE1=zeros(1,Lch); 
       hE1(Delay+1)=channel; % cir: channel impulse response
       hEr1=hE1;
       HEr1=fft([hEr1 zeros(1,Nfft-Lch)]); % Channel frequency response
      

       %% Round 2 -Antenna 2
        %h_ab  a=usr b=antena
   
      %user1  
       channel=(randn(1,Ntap)+j*randn(1,Ntap)).*sqrt(Power/2);%%%%%
       h12=zeros(1,Lch); 
       h12(Delay+1)=channel; % cir: channel impulse response
       h12r2=h12;
       H12r2=fft([h12r2 zeros(1,Nfft-Lch)]); % Channel frequency response
      
       %user2
      % rng;
       channel=(randn(1,Ntap)+j*randn(1,Ntap)).*sqrt(Power/2);%%%%%
       h22=zeros(1,Lch); 
       h22(Delay+1)=channel; % cir: channel impulse response
       h22r2=h22;
       H22r2=fft([h22r2 zeros(1,Nfft-Lch)]); % Channel frequency response
     
       %Eve
       channel=(randn(1,Ntap)+j*randn(1,Ntap)).*sqrt(Power/2);%%%%%
       hE2=zeros(1,Lch); 
       hE2(Delay+1)=channel; % cir: channel impulse response
       hEr2=hE2;
       HEr2=fft([hEr2 zeros(1,Nfft-Lch)]); % Channel frequency response
      

       %%%
      % User1 Tx______________________________________________________________
      X= randi([0,M-1], Nused*Nframe,1)'; % bit: integer vector
      %   X= randint(1,Nused*Nframe,M); % bit: integer vector %Non binary data 0 1 2 3
      Xmod= qammod(X,M,'gray')/norms(Nbps); %%
      if NgType~=2, x_GI=zeros(1,Nframe*Nsym);
      elseif NgType==2, x_GI= zeros(1,Nframe*Nsym+Ng);
        % Extend an OFDM symbol by Ng zeros 
      end
      
      % User2 Tx______________________________________________________________
    
      X2= randi([0,M-1], Nused*Nframe,1)'; % bit: integer vector
      %   X= randint(1,Nused*Nframe,M); % bit: integer vector %Non binary data 0 1 2 3
      Xmod2= qammod(X2,M,'gray')/norms(Nbps); %%
      if NgType~=2, x_GI2=zeros(1,Nframe*Nsym);
      elseif NgType==2, x_GI2= zeros(1,Nframe*Nsym+Ng);
        % Extend an OFDM symbol by Ng zeros 
      end

%%Modulted signal for user1 and user2      
      U1=Xmod;
      U2=Xmod2;
 
  %Precoding matrix 
  %h_ab  a=usr b=antena
   
%  M_1a=diag(H22r2)*(diag(H22r2)*diag(H11r1)-diag(H21r1)*diag(H12r2))^-1;
%  M_1b=-(M_1a*diag(H21r1))/diag(H22r2);
% 
%  M_2a=diag(H12r2)*(diag(H21r1)*diag(H12r2)-diag(H22r2)*diag(H11r1))^-1;
%  M_2b=-(M_2a*diag(H11r1))/diag(H12r2);
%   
%  
 M_1a=diag(H22r2)*(diag(H22r2)*diag(H11r1)-diag(H21r1)*diag(H12r2))^-1;
 M_1b=-(M_1a*diag(H21r1))/diag(H22r2);

 M_2a=diag(H12r2)*(diag(H21r1)*diag(H12r2)-diag(H22r2)*diag(H11r1))^-1;
 M_2b=-(M_2a*diag(H11r1))/diag(H12r2);
 
 
%  % Precoder matrices
%        A = diag(H11r1);
%        B = diag(H12r2);
%        C = diag(H21r1);
%        D = diag(H22r2);
%        
%        p1=-D/(C*B-D*A);
%        p3=C/(C*B-D*A);
% 
%        p2=-B/(D*A-C*B);
%        p4=A/(D*A-C*B);
%        
%   %%    
% %       Transmiters
%       Tx1=(1/sqrt(2))*X1mod*(p1.*1/sqrt(abs(p1))) + (1/sqrt(2))*X2mod*(p2.*1/sqrt(abs(p2)));
%       %Antenna 2    
%       Tx2=(1/sqrt(2))*X1mod*(p3.*1/sqrt(abs(p3))) + (1/sqrt(2))*X2mod*(p4.*1/sqrt(abs(p4)));
%       p1nom = (norm(diag(p1),2))^2;
%       p2nom = (norm(diag(p2),2))^2;
%       p3nom = (norm(diag(p3),2))^2;
%       p4nom = (norm(diag(p4),2))^2;
%       
%       Tx1=(1/sqrt(2))*U1*(p1*1/sqrt(p1nom)) + (1/sqrt(2))*U2*(p2*1/sqrt(p2nom));
%       %Antenna 2    
%       Tx2=(1/sqrt(2))*U1*(p3*1/sqrt(p3nom)) + (1/sqrt(2))*U2*(p4*1/sqrt(p4nom));
  

 %% Normaization factore
     snr = EbN0(i)+10*log10(Nbps*(Nused/Nfft)); % SNR vs. Eb/N0
     noise_var = Nfft*0.5*10^(-snr/10); 
  
    % EbN0
 %%Super position
Tx1=1/sqrt(2)*U1*M_1a+1/sqrt(2)*U2*M_2a;
 %% Antenna 2    
Tx2=1/sqrt(2)*U1*M_1b+1/sqrt(2)*U2*M_2b;

 %%Super position
% Tx1=sqrt(1)*U1*M_1a+sqrt(1)*U2*M_2a;
%  %% Antenna 2    
% Tx2=sqrt(1)*U1*M_1b+sqrt(1)*U2*M_2b;

     x1np= ifft(U1);
     x_GI1np= guard_interval(Ng,Nfft,NgType,x1np);       
          
         %Non precoded u2
     x2np= ifft(U2);
     x_GI2np= guard_interval(Ng,Nfft,NgType,x2np);       
 
%% PRecoded: This is imp 
         x1= ifft(Tx1)/norm(ifft(Tx1),2);  %% This is how i am normalizign it :step1
         norm(ifft(Tx1),2);
        % pause()
         x_GI1= guard_interval(Ng,Nfft,NgType,x1);       
         y11= conv(x_GI1,h11r1); %u1 r1
         y21=(conv(x_GI1,h21r1)); %u2 r1
         yE1= conv(x_GI1,hEr1); %u2 r1
         %EbN0=20;
         snr = EbN0(i)+10*log10(Nbps*(Nused/Nfft)); % SNR vs. Eb/N0
        %snr=snr/2;
         EEs=sum(x1).^2/length(x1);    
         %EEs=1;
         noise_mag  = sqrt(0.5*EEs*((10.^(-snr/10))));
         y_G11r1 = norm(ifft(Tx1),2)*y11+ noise_mag*(randn(size(y11))+j*randn(size(y11))); %% This is how i am normalizign : step2
         y_G21r1 = norm(ifft(Tx1),2)*y21+ noise_mag*(randn(size(y21))+j*randn(size(y21)));
         y_G21E1 = norm(ifft(Tx1),2)*yE1+ noise_mag*(randn(size(yE1))+j*randn(size(yE1)));
 
         
         Y11= fft(remove_GI(Ng,Nsym,NgType,y_G11r1));
         Y21= fft(remove_GI(Ng,Nsym,NgType,y_G21r1));
         YE1= fft(remove_GI(Ng,Nsym,NgType,y_G21E1));
 
         %% %
 %A2   
    x2= ifft(Tx2)/norm(ifft(Tx2),2);
    x_GI2= guard_interval(Ng,Nfft,NgType,x2);       
    y12=conv(x_GI2,h12r2); %u1 r2
    y22= conv(x_GI2,h22r2); %u2 r2
    yE2= conv(x_GI2,hEr2); %u2 r1    
 
    
     %% PAPR
         
         if i==1
         PAPRtx1(1,m)=max((abs(x_GI1)).^2)/mean((abs(x_GI1)).^2);
         PAPRtx2(1,m)=max((abs(x_GI2)).^2)/mean((abs(x_GI2)).^2);
         PAPRtx3(1,m)=max((abs(x_GI1np)).^2)/mean((abs(x_GI1np)).^2);
         PAPRtx4(1,m)=max((abs(x_GI2np)).^2)/mean((abs(x_GI2np)).^2);
         %pause();
          end
      
    
         %EbN0=10;
  
       EEs=sum(x2).^2/length(x2);    
       noise_mag  = sqrt(0.5*EEs*((10.^(-snr/10))));
    % noise_mag  = sqrt(0.5*EEs*10.^(-snr/10));
         
       y_G12r2 = norm(ifft(Tx2),2)* y12+ noise_mag*(randn(size(y12))+j*randn(size(y12)));
       y_G22r2 = norm(ifft(Tx2),2)* y22+ noise_mag*(randn(size(y22))+j*randn(size(y22)));
       y_G21E2 = norm(ifft(Tx2),2)*yE2+ noise_mag*(randn(size(yE2))+j*randn(size(yE2)));
 
       Y12= fft(remove_GI(Ng,Nsym,NgType,y_G12r2));
       Y22= fft(remove_GI(Ng,Nsym,NgType,y_G22r2));
       YE2= fft(remove_GI(Ng,Nsym,NgType,y_G21E2));
         
         %%
         YU1= (Y11)+(Y12);
         YU2= (Y21)+(Y22);
         YUe= YE1+YE2;
%          figure(5);stem(abs(YU1)); hold on; stem(abs(U1));%stem(Tx1)
%          figure(6);stem(abs(YU2)); hold on; stem(abs(U2));%stem(Tx1)
%          figure(7);stem(abs(YUe)); hold on; stem(abs(U1));%stem(Tx1)
%          figure(8);stem(abs(YUe)); hold on; stem(abs(U2));%stem(Tx1)

        UU1i=qamdemod(YU1*norms(Nbps),M,'gray');
        UU2i=qamdemod(YU2*norms(Nbps),M,'gray');
  
        UUee=qamdemod(YUe*norms(Nbps),M,'gray');
  
    
%figure(31);stem(UU1i); hold on; stem(X);%stem(Tx1)
user1n=length(find(UU1i==X));

%figure(32);stem(UU2i); hold on; stem(X2);%stem(Tx1)
user2n=length(find(UU2i==X2));

%% External Eve
user1en=length(find(UUee==X));
user2en=length(find(UUee==X2));
%% Internal Eve
%%%
userI1en=length(find(X2==X));
userI2en=length(find(X==X2));


SERu1=SERu1+(Nfft-user1n);
SERu2=SERu2+(Nfft-user2n);
SERu1E=SERu1E+(Nfft-user1en);
SERu2E=SERu2E+(Nfft-user2en);
%% I
SERuI1E=SERuI1E+(Nfft-userI1en);
SERuI2E=SERuI2E+(Nfft-userI2en);


symtot=symtot+Nfft;
%%
 BERu1=BERu1+sum(sum(de2bi(UU1i,Nbps)~=de2bi(X,Nbps)));
      BERu1p(i,m)=sum(sum(de2bi(UU1i,Nbps)~=de2bi(X,Nbps))); %for PER
 BERu2= BERu2+sum(sum(de2bi(UU2i,Nbps)~=de2bi(X2,Nbps)));
      BERu2p(i,m)=sum(sum(de2bi(UU2i,Nbps)~=de2bi(X2,Nbps)));

      %% Throughput
     TERu1=TERu1+(Nfft*Nbps)-sum(sum(de2bi(UU1i,Nbps)~=de2bi(X,Nbps)));     
     TERu2=TERu2+(Nfft*Nbps)-sum(sum(de2bi(UU2i,Nbps)~=de2bi(X2,Nbps)));
   TERu1E=TERu1E+(Nfft*Nbps)-sum(sum(de2bi(UUee,Nbps)~=de2bi(X,Nbps)));
   TERu2E=TERu2E+(Nfft*Nbps)-sum(sum(de2bi(UUee,Nbps)~=de2bi(X2,Nbps)));
   TERI1E=TERI1E+(Nfft*Nbps)-sum(sum(de2bi(X,Nbps)~=de2bi(X2,Nbps)));
   TERI2E=TERI2E+(Nfft*Nbps)-sum(sum(de2bi(X2,Nbps)~=de2bi(X,Nbps)));
      %% External Eve
  BERu1E=BERu1E+sum(sum(de2bi(UUee,Nbps)~=de2bi(X,Nbps)));
      BERu1Ep(i,m)=sum(sum(de2bi(UUee,Nbps)~=de2bi(X,Nbps)));
 BERu2E=BERu2E+sum(sum(de2bi(UUee,Nbps)~=de2bi(X2,Nbps)));
     BERu2Ep(i,m)=sum(sum(de2bi(UUee,Nbps)~=de2bi(X2,Nbps)));
 
     %% Internal Eve
 BERI1E=BERI1E+sum(sum(de2bi(X,Nbps)~=de2bi(X2,Nbps)));
  BERI1Ep(i,m)=sum(sum(de2bi(X,Nbps)~=de2bi(X2,Nbps)));
 BERI2E=BERI2E+sum(sum(de2bi(X2,Nbps)~=de2bi(X,Nbps)));
  BERI2Ep(i,m)=sum(sum(de2bi(X2,Nbps)~=de2bi(X,Nbps)));

 bittot=bittot+Nfft*Nbps; 

%[user1n user2n user1en user2en]

   end
   
   
%SER
SERu11(i)=SERu1/symtot;
SERu21(i)=SERu2/symtot;
SERu1E1(i)=SERu1E/symtot;
SERu2E1(i)=SERu2E/symtot;
%symtot1=symtot+Nfft;
SERuI1Ei(i)=SERuI1E/symtot;
SERuI2Ei(i)=SERuI2E/symtot;
 
%BER
 BERu11(i)=BERu1/bittot;
 BERu21(i)=BERu2/bittot;
 BERu1E1(i)=BERu1E/bittot;
 BERu2E1(i)=BERu2E/bittot;

 BERI1Ei(i)=BERI1E/bittot;
 BERI2Ei(i)=BERI2E/bittot;

 % Thoruput
 resour=N_iter*Nfft;
 TERu11(i)=TERu1/(resour);
 TERu21(i)=TERu2/resour;
 TERu1E1(i)=TERu1E/resour;
 TERu2E1(i)=TERu2E/resour;
 TERI1Ei(i)=TERI1E/resour;
 TERI2Ei(i)=TERI2E/resour;
 
 % Thoruput : cosndeing half resouces
 resour=N_iter*Nfft;
 TERu11h(i)=TERu1/(resour*0.5);
 TERu21h(i)=TERu2/(resour*0.5);
 TERu1E1h(i)=TERu1E/(resour*0.5);
 TERu2E1h(i)=TERu2E/(resour*0.5);
 TERI1Eih(i)=TERI1E/(resour*0.5);
 TERI2Eih(i)=TERI2E/(resour*0.5);
 
 
 %PER
 PERu11(i,:)=BERu1p(i,:)>0;
 PERu11p(i)=sum(PERu11(i,:)')/N_iter;
 
 PBERu21(i,:)=BERu2p(i,:)>0;
 PBERu21p(i)=sum(PBERu21(i,:)')/N_iter;
 
 PBERu1E1(i,:)=BERu1Ep(i,:)>0;
 PBERu1E1p(i)=sum(PBERu1E1(i,:)')/N_iter;
 
 PBERu2E1(i,:)=BERu2Ep(i,:)>0;
 PBERu2E1p(i)=sum(PBERu2E1(i,:)')/N_iter;
 
 %%Internal
 PBERI1E1(i,:)=BERI1Ep(i,:)>0;
 PBERI1E1p(i)=sum(PBERI1E1(i,:)')/N_iter;
 
 PBERI2E1(i,:)=BERI2Ep(i,:)>0;
 PBERI2E1p(i)=sum(PBERI2E1(i,:)')/N_iter;
 
% PERu1=
% PERu2=
% PERu1E=
% PERu2E=

   
end
%%
 [BERu11DSISO, BERu11SISO, TERu11DSISO, TERu1_SISO] = WP2_BER_DSISO(N_iter);
%%
figure(1);
semilogy(EbN0,TERu11,'or-','LineWidth',2); hold on;
semilogy(EbN0,TERu21,'*b-','LineWidth',2);hold on;
% semilogy(EbN0,TERu11DSISO,'vc-','LineWidth',2); hold on;
semilogy(EbN0,TERu1_SISO,'mp-','LineWidth',2); 
grid on;
xlabel('SNR (dB)');
ylabel('Throughput rate');
legend('RX1-proposed','RX2-proposed','Conventional-2 by 2 MIMO','fontweight','bold');

% PAPRtx1dB=10*log10(PAPRtx1); % in dB
% [F1,X1] = ecdf(PAPRtx1dB); %CDF
% 
% PAPRtx2dB=10*log10(PAPRtx2); % in dB
% [F2,X2] = ecdf(PAPRtx2dB); %CDF
% % 
% end

