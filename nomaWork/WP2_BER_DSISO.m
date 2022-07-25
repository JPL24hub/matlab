
function[BERu11DSISO, BERu11SISO, TERu11DSISO, TERu1_SISO] = WP2_BER_DSISO(N_iter)
%%
warning off;

Nbps=1; 
M=2^Nbps;           % Modulation order=2/4/6 for QPSK/16QAM/64QAM
EbN0=1:1:30;   % EbN0
% N_iter=1000;        % no of packets   % Number of iterations for each EbN0
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
     BERu1DSISO=0;
     BERu2DSISO=0;
     BERu1SISO=0;
     
     TERu1DSISO=0;
     TERu2DSISO=0;
     TERu1SISO=0;
     
     bittot=0;
     
     for m=1:N_iter
       %% antenna 1 USER 1
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
       %% Generate Signal for each user/tx, modulate, and add CP
       
      % User1 Tx______________________________________________________________
      X= randi([0,M-1], Nused*Nframe,1)'; % bit: integer vector (1)
      Xmod= qammod(X,M,'gray')/norms(Nbps); 	
       % User2 Tx______________________________________________________________
%       X2= randi([0,M-1], Nused*Nframe,1)'; % bit: integer vector
%       Xmod2= qammod(X2,M,'gray')/norms(Nbps);
      %%Modulted signal for user1 and user2      
      U1=Xmod;
%       U2=Xmod2;
       
       %% Super position
        %%Antenna 1
       Tx1=sqrt(1)*U1;
       Tx2=sqrt(1)*U1;
           
       %% PRecoded: This is imp (from X to Y)
       %Tx1 ROUND 1
         x1= ifft(Tx1)/norm(ifft(Tx1),2);  %% This is how i am normalizign it :step1
         norm(ifft(Tx1),2);
         
        % pause()
         x_GI1= guard_interval(Ng,Nfft,NgType,x1);       
         y11= conv(x_GI1,h11r1); %u1 r1
         y21=(conv(x_GI1,h21r1)); %u2 r1
         
         %EbN0=20;
         snr = EbN0(k)+10*log10(Nbps*(Nused/Nfft)); % SNR vs. Eb/N0
         EEs=sum(x1).^2/length(x1);    
         %EEs=1;
         noise_mag  = sqrt(0.5*EEs*((10.^(-snr/10))));
         y_G11r1 = norm(ifft(Tx1),2)*y11+ noise_mag*(randn(size(y11))+1i*randn(size(y11))); %% This is how i am normalizign : step2
         y_G21r1 = norm(ifft(Tx1),2)*y21+ noise_mag*(randn(size(y21))+1i*randn(size(y21)));
         
         Y11= fft(remove_GI(Ng,Nsym,NgType,y_G11r1));
         Y21= fft(remove_GI(Ng,Nsym,NgType,y_G21r1));

         %% Tx2  ROUND 2
         
         x2= ifft(Tx2)/norm(ifft(Tx2),2);
         x_GI2= guard_interval(Ng,Nfft,NgType,x2);       
         y12=conv(x_GI2,h12r2); %u1 r2
         y22= conv(x_GI2,h22r2); %u2 r2
         
         EEs=sum(x2).^2/length(x2);    
        noise_mag  = sqrt(0.5*EEs*((10.^(-snr/10))));  
        y_G12r2 = norm(ifft(Tx2),2)* y12+ noise_mag*(randn(size(y12))+1i*randn(size(y12)));
        y_G22r2 = norm(ifft(Tx2),2)* y22+ noise_mag*(randn(size(y22))+1i*randn(size(y22)));
        
        Y12= fft(remove_GI(Ng,Nsym,NgType,y_G12r2));
        Y22= fft(remove_GI(Ng,Nsym,NgType,y_G22r2));
         
         %% MRC
         YU1= (H11r1)'.*transpose(Y11)+(H12r2)'.*transpose(Y12);
         YU2= (H21r1)'.*transpose(Y21)+(H22r2)'.*transpose(Y22);
         YU1SISO= (H11r1)'.*transpose(Y11);
         
         %%Demodulate
         UU1i=qamdemod(YU1*norms(Nbps),M,'gray');
         UU2i=qamdemod(YU2*norms(Nbps),M,'gray');
         UU1iSISO=qamdemod(YU1SISO*norms(Nbps),M,'gray');
         

          %% BER
          
          BERu1DSISO= BERu1DSISO+sum(sum(de2bi(UU1i,Nbps)~=de2bi(X,Nbps)));
          BERu2DSISO= BERu2DSISO+sum(sum(de2bi(UU2i,Nbps)~=de2bi(X,Nbps)));
          BERu1SISO= BERu1SISO+sum(sum(de2bi(UU1iSISO,Nbps)~=de2bi(X,Nbps)));
       
          TERu1DSISO=TERu1DSISO+(Nfft*Nbps)-sum(sum(de2bi(UU1i,Nbps)~=de2bi(X,Nbps))); 
          TERu2DSISO=TERu2DSISO+(Nfft*Nbps)-sum(sum(de2bi(UU2i,Nbps)~=de2bi(X,Nbps)));
          TERu1SISO=TERu1SISO+(Nfft*Nbps)-sum(sum(de2bi(UU1iSISO,Nbps)~=de2bi(X,Nbps)));
          
          
          bittot=bittot+Nfft*Nbps; 
         
     end
     
     %BER
      BERu11DSISO(k)=BERu1DSISO/bittot;
      BERu21DSISO(k)=BERu2DSISO/bittot;
      BERu11SISO(k)=BERu1SISO/bittot;
      
      resour=N_iter*Nfft;
      TERu11DSISO(k)=TERu1DSISO/(resour);
      TERu21DSISO(k)=TERu2DSISO/resour;
      TERu1_SISO(k)=TERu1SISO/resour;

      
end
%%
end
%BER
% figure(1);
% semilogy(EbN0,BERu11DSISO,'*-','LineWidth',1); hold on;
% semilogy(EbN0,BERu11SISO,'s-','LineWidth',1);
% xlabel('SNR (dB)');
% ylabel('Bit Error Rate(BER)');
% legend('BER (SISO-dual)','BER (SISO)');
% %%
% figure(3);
% semilogy(EbN0,TERu11DSISO,'s-','LineWidth',1);hold on;
% % semilogy(EbN0,TERu21DSISO,'*b--','LineWidth',2);hold on;
% semilogy(EbN0,TERu1_SISO,'*-','LineWidth',1);hold on;
% xlabel('SNR (dB)');
% ylabel('Throughput (BER)');
% legend('Throughput (SISO)','Throughput (SISO-dual)');