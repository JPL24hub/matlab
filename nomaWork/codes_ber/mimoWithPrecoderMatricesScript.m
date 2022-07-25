%New Wireless Multiple Input Multiple Output Communication Technique with Enhance Security Using Precoder Matrices
warning off;
clear all
clc


Nbps=1; 
M=2^Nbps;           % Modulation order=2/4/6 for QPSK/16QAM/64QAM
EbN0=0:1:40;        % EbN0, SNR
%%
d=0:1:80;
di=length(d);
N_iter=6000;        % no of packets   % Number of iterations for each EbN0
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

NgType=1;          % NgType=1/2 for cyclic prefix/zero padding
if NgType==1, nt='CP';  elseif NgType==2, nt='ZP';   end

 PERu11=[];
 PERu11p=[];
      
 PERu21=[];
 PERu21p=[];
      
 PBERu1E1=[];
 PBERu1E1p=[];
      
PBERu2E1=[];
PBERu2E1p=[];

for k=1:length(EbN0)
     BER_Rx1=0;
     BER_Rx2=0;
     BER_EVE_Rx1=0;
     BER_EVE_Rx2=0;
     BER_all_antennas=0;  
     
     TER_Rx1=0;
     TER_Rx2=0;
     TER_EVE_Rx1=0;
     TER_EVE_Rx2=0;
   
     bittot=0;
     newbittot=0;
     NR=2;
     
     for m=1:N_iter
      %h_ab  a=Tx b=Rx
% Tx1 to to all Rx
      %Tx1 to Rx1 
      channel11=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);
       h11=zeros(1,Lch); 
       h11(Delay+1)=channel11; % cir: channel impulse response
       h11Tx1=h11;
       H11r1=fft([h11Tx1 zeros(1,Nfft-Lch)]); % Channel frequency response
       
      %Tx1 to Rx2
       channel21=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);
       h21=zeros(1,Lch); 
       h21(Delay+1)=channel21; % cir: channel impulse response
       h21Tx1=h21;
       H21Tx1=fft([h21Tx1 zeros(1,Nfft-Lch)]); % Channel frequency response
      
       %Tx1 to EVE
       channelE1= channel11 + (randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);%*d(di);
       hE1=zeros(1,Lch); 
       hE1(Delay+1)=channelE1; % cir: channel impulse response
       hETx1=hE1;
       HEr1=fft([hETx1 zeros(1,Nfft-Lch)]); % Channel frequency response
       
% Tx2 to all Rx

      %Tx2 to Rx1
       channel12=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);	
       h12=zeros(1,Lch); 
       h12(Delay+1)=channel12; % cir: channel impulse response
       h12Tx2=h12;
       H12r2=fft([h12Tx2 zeros(1,Nfft-Lch)]); % Channel frequency response
      
      %Tx2 to Rx2
       channel22=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);	
       h22=zeros(1,Lch); 
       h22(Delay+1)=channel22; % cir: channel impulse response
       h22Tx2=h22;
       H22r2=fft([h22Tx2 zeros(1,Nfft-Lch)]); % Channel frequency response
     
       %Tx2 to EVE
       channelE2= channel12 + (randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);%*d(di);	
       hE2=zeros(1,Lch); 
       hE2(Delay+1)=channelE2; % cir: channel impulse response
       hETx2=hE2;
       HEr2=fft([hETx2 zeros(1,Nfft-Lch)]); % Channel frequency response
       
%  Generate Signal for each antenna, modulate, and add CP
      % X1 Data
      X1= randi([0,M-1], Nused*Nframe,1)'; % bit: integer vector (1)
      Xmod= qammod(X1,M,'gray')/norms(Nbps); 	
      
      % X2 Data
      X2= randi([0,M-1], Nused*Nframe,1)'; % bit: integer vector
      Xmod2= qammod(X2,M,'gray')/norms(Nbps);
           
      U1=Xmod;
      U2=Xmod2;    
      
% Precoder matrices
       A = diag(H11r1);
       B = diag(H12r2);
       C = diag(H21Tx1);
       D = diag(H22r2);
       
       p1=-D/(C*B-D*A);
       p3=C/(C*B-D*A);

       p2=-B/(D*A-C*B);
       p4=A/(D*A-C*B);
       
% Superposition
%       %Antenna 1
%       Tx1=sqrt(2)*U1*(p1)+sqrt(2)*U2*(p2);
%       %Antenna 2    
%       Tx2=sqrt(2)*U1*(p3)+sqrt(2)*U2*(p4);
      
            %Antenna 1
%       Tx1=(1/sqrt(2))*U1*(p1.*1/norm(diag(p1),2)) + (1/sqrt(2))*U2*(p2.*1/norm(diag(p2),2));
%       %Antenna 2    
%       Tx2=(1/sqrt(2))*U1*(p3.*1/norm(diag(p3),2)) + (1/sqrt(2))*U2*(p4.*1/norm(diag(p4),2));
      
      Tx1=(1/sqrt(2))*U1*(p1.*1/sqrt(abs(p1))) + (1/sqrt(2))*U2*(p2.*1/sqrt(abs(p2)));
      %Antenna 2    
      Tx2=(1/sqrt(2))*U1*(p3.*1/sqrt(abs(p3))) + (1/sqrt(2))*U2*(p4.*1/sqrt(abs(p4)));
      
      %% FOR PAPR    
      x1np= ifft(U1);
      x_GI1np= guard_interval(Ng,Nfft,NgType,x1np);       
      x2np= ifft(U2);
      x_GI2np= guard_interval(Ng,Nfft,NgType,x2np);
      
      np1=sqrt(1/(mean(trace(p1*p1'))));
      np2=sqrt(1/(mean(trace(p2*p2'))));
      np3=sqrt(1/(mean(trace(p3*p3'))));
      np4=sqrt(1/(mean(trace(p4*p4'))));
      
      
% Antenna 1 (Tx1)
%         x1= ifft(Tx1)*(np1+np2);  %% This(np3+np4) is how i am normalizign it :step1
%         x1=ifft(Tx1)/norm(ifft(Tx1),2);
        x1=ifft(Tx1);
        x1WithGI=guard_interval(Ng,Nfft,NgType,x1);       
        y11=conv(x1WithGI,h11Tx1); %u1 r1
        y21=conv(x1WithGI,h21Tx1); %u1 r1
        y31=conv(x1WithGI,hETx1); %u1 r1
        
        
        % Antenna 2 (Tx2)
%         x2= ifft(Tx2)*(np3+np4);
%         x2= ifft(Tx2)/norm(ifft(Tx2),2);
        x2=ifft(Tx2);
        x2WithGI=guard_interval(Ng,Nfft,NgType,x2);
        
        fadedSig=h_t*txSig;
        
        y12=conv(x2WithGI,h12Tx2); %u1 r2
        y22=conv(x2WithGI,h22Tx2); %u1 r2
        y32=conv(x2WithGI,hETx2);  %u1 r2
        
         
         %EbN0=20;
        %Noise 1
        snr = EbN0(k)+10*log10(Nbps*(Nused/Nfft)); % SNR vs. Eb/N0         
        EEs=(sum(x1).^2)/length(x1);
%         EEs=1;
        noise_mag  = sqrt(0.5*EEs*((10.^(-snr/10))));
        
        n1 =  noise_mag*(randn(size(y11))+1i*randn(size(y11)));

       %Noise 2
       EEs=(sum(x2).^2)/length(x2);    
       noise_mag = sqrt(0.5*EEs*((10.^(-snr/10))));    
       n2 = noise_mag*(randn(size(y12))+1i*randn(size(y12)));
        
        %Noise 3
       EEs=(sum(x1).^2)/length(x2);    
       noise_mag = sqrt(0.5*EEs*((10.^(-snr/10))));    
       n3 = noise_mag*(randn(size(y32))+1i*randn(size(y32)));
       
        % %%%%%%% FOR PAPR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         if k==1
         PAPRtx1(1,m)=max((abs(x1WithGI)).^2)/mean((abs(x1WithGI)).^2);
         PAPRtx2(1,m)=max((abs(x2WithGI)).^2)/mean((abs(x2WithGI)).^2);
         PAPRtx3(1,m)=max((abs(x_GI1np)).^2)/mean((abs(x_GI1np)).^2);
         PAPRtx4(1,m)=max((abs(x_GI2np)).^2)/mean((abs(x_GI2np)).^2);
         %pause();
         end
         
         % %%%%%%% FOR PAPR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    %==============Received signals=================================
%        Y1 = norm(ifft(Tx1),2)*y11 + norm(ifft(Tx2),2)*y12 + n1; 
         Y1 = y11 + y12 + n1;       
%        Y2 = norm(ifft(Tx1),2)*y21 + norm(ifft(Tx2),2)*y22 + n2;
         Y2 = y21 + y22 + n2;    
       Y3 = norm(ifft(Tx1),2)*y31 + norm(ifft(Tx2),2)*y32 + n3;
       
       Y1new = fft(remove_GI(Ng,Nsym,NgType,Y1));
       Y2new = fft(remove_GI(Ng,Nsym,NgType,Y2));
       Y3new = fft(remove_GI(Ng,Nsym,NgType,Y3));
        
        
% Demodulation
        Rx1=qamdemod(Y1new*norms(Nbps),M,'gray');
        Rx2=qamdemod(Y2new*norms(Nbps),M,'gray');
        RxEVE=qamdemod(Y3new*norms(Nbps),M,'gray');

        %% BER 
        BER_Rx1= BER_Rx1+sum(sum(de2bi(Rx1,Nbps)~=de2bi(X1,Nbps)));
        BER_Rx2= BER_Rx2+sum(sum(de2bi(Rx2,Nbps)~=de2bi(X2,Nbps)));
        BER_EVE_Rx1= BER_EVE_Rx1+sum(sum(de2bi(RxEVE,Nbps)~=de2bi(X1,Nbps)));
        BER_EVE_Rx2= BER_EVE_Rx2+sum(sum(de2bi(RxEVE,Nbps)~=de2bi(X2,Nbps)));
        
        BER_all_antennas = BER_all_antennas+sum(sum(de2bi([Rx1;Rx2],Nbps)~=de2bi([X1';X2'],Nbps)));
        
        %% Throughput
        TER_Rx1=TER_Rx1+(Nfft*Nbps)-sum(sum(de2bi(Rx1,Nbps)~=de2bi(X1,Nbps))); 
        TER_Rx2=TER_Rx2+(Nfft*Nbps)-sum(sum(de2bi(Rx2,Nbps)~=de2bi(X2,Nbps)));
        TER_EVE_Rx1=TER_EVE_Rx1+(Nfft*Nbps)-sum(sum(de2bi(RxEVE,Nbps)~=de2bi(X1,Nbps)));
        TER_EVE_Rx2=TER_EVE_Rx2+(Nfft*Nbps)-sum(sum(de2bi(RxEVE,Nbps)~=de2bi(X2,Nbps)));
        
        %% PAPR
        BERu1p(k,m)=sum(sum(de2bi(Rx1,Nbps)~=de2bi(X1,Nbps))); %for PER
        BERu2p(k,m)=sum(sum(de2bi(Rx2,Nbps)~=de2bi(X2,Nbps))); %for PER
        BERu1Ep(k,m)=sum(sum(de2bi(RxEVE,Nbps)~=de2bi(X1,Nbps)));
        BERu2Ep(k,m)=sum(sum(de2bi(RxEVE,Nbps)~=de2bi(X2,Nbps)));


        
        bittot=bittot+Nfft*Nbps;
        newbittot=newbittot+2*Nfft*Nbps;
        
 
        
        c(m)=(NR*log(1+(10^(snr/10))*abs(normrnd(0,1)))/log(2));
        
    
     end
     di=di-1;
     
    %CAPACITY 
     SNR(k)=snr;
     ysnr(k)=mean(c);
     shanonCapacity(k)=(log(1+10^(snr/10)))/log(2);
     
     
     %BER
      BER(k) = BER_all_antennas/newbittot;
      BER_Tx1Rx1(k)=BER_Rx1/bittot;         % BER for TX1 TO RX1
      BER_Tx2Rx2(k)=BER_Rx2/bittot;         % BER for TX2 TO RX2
      BER_EVE_Tx1(k)=BER_EVE_Rx1/bittot;    % BER for EVE TO RX1
      BER_EVE_Tx2(k)=BER_EVE_Rx2/bittot;    % BER for EVE TO RX2
      
      %Thoruput
      resour=N_iter*Nfft;
      TER_Tx1Rx1(k)=TER_Rx1/(resour);
      TER_Tx2Rx2(k)=TER_Rx2/resour;
      TER_EVE_Tx1(k)=TER_EVE_Rx1/resour;
      TER_EVE_Tx2(k)=TER_EVE_Rx2/resour;
      
      %PER
      PERu11(k,:)=BERu1p(k,:)>0;
      PERu11p(k)=sum(PERu11(k,:)')/N_iter;
      
      PERu21(k,:)=BERu2p(k,:)>0;
      PERu21p(k)=sum(PERu21(k,:)')/N_iter;%Creating
      
      PBERu1E1(k,:)=BERu1Ep(k,:)>0;
      PBERu1E1p(k)=sum(PBERu1E1(k,:)')/N_iter;
      
      PBERu2E1(k,:)=BERu2Ep(k,:)>0;
      PBERu2E1p(k)=sum(PBERu2E1(k,:)')/N_iter;
end
%%
% [Eb_N0_dB,theoryBer_nRx1,theoryBerMRC_nRx2,simBer]=MIMOwithZF();
% [EbN01,BER_simm]=OFDMberINawgn();
% [EbN02,BER_simimp]=OFDMberINawgnImpafect();


%% BER
semilogy(EbN01,BER_simm,'om-','LineWidth',2); hold on
semilogy(EbN01,BER_simm,'*g--','LineWidth',2); hold on
semilogy(EbN0,BER_Tx1Rx1,'ok--','LineWidth',2); hold on
semilogy(EbN0,BER_Tx2Rx2,'*r--','LineWidth',2); hold on
semilogy(EbN0,BER_Tx1Rx1_DIV,'ok-','LineWidth',2); hold on
semilogy(EbN0,BER_Tx2Rx2_DIV,'*m--','LineWidth',2); hold on
semilogy(Eb_N0_dB,theoryBer_nRx1,'pb-','LineWidth',2); hold on
semilogy(Eb_N0_dB,simBer,'pc--','LineWidth',2); hold on

ax = gca;
xlabel('SNR (dB)');
ax.XAxis.LineWidth = 2;
ylabel('Bit Error Rate(BER)');
ax.YAxis.LineWidth = 2;
grid on;
legend('Rx1 proposed-sim', 'Rx2 proposed-sim','Rx1-proposed','Rx2-proposed','Rx1-with-Dvsty-proposed','Rx2-with-Dvsty-proposed','theory (nTx=1,nRx=1)','sim (nTx=2, nRx=2, ZF)','fontweight','bold');
%%




































% Plots
%% CAPACITY
% 
% plot(SNR, ysnr,'--pm','LineWidth',2);hold on;
% plot(SNR, shanonCapacity, '*r--', 'LineWidth',2); hold off;
% ylabel('Capacity (bit/s/Hz)');
% xlabel('SNR (dB)');
% legend('Proposed Model', 'Shannon','fontweight','bold');
% grid on;



% 
% %%
% 
% % figure(1);
% 
% semilogy(EbN0(1,1:col),BER_Tx1Rx1(1,1:col),'ok-','LineWidth',2); hold on
% semilogy(EbN0(1,1:col),BER_Tx2Rx2(1,1:col),'*r--','LineWidth',2); hold on
% semilogy(EbN0(1,1:col),BER_EVE_Tx1(1,1:col),'*g-','LineWidth',2);hold on
% semilogy(EbN0(1,1:col),BER_EVE_Tx2(1,1:col),'ob--','LineWidth',2); hold on;
% %%
% semilogy(EbN0(1,1:col),BER_simimp(1,1:col),'-pm','LineWidth',2); hold on; % RUN MIMOwihtDIV
% semilogy(EbN0(1,1:col),BER_simimp(1,1:col),'--pg','LineWidth',2); hold on; % RUN MIMOwihtDIV
% semilogy(EbN0(1,1:col),BER_EVE_Tx1(1,1:col),'-pk','LineWidth',2); hold on; % RUN MIMOwihtDIV
% semilogy(EbN0(1,1:col),BER_EVE_Tx2(1,1:col),'--pr','LineWidth',2); hold on; % RUN MIMOwihtDIV
% semilogy(EbN0(1,1:col),theoryBer_nRx1(1,1:col),'-pb','LineWidth',2); hold on;
% semilogy(EbN0(1,1:col),simBer(1,1:col),'--pc','LineWidth',2); hold on;
% % semilogy(EbN0(1,1:col),BER(1,1:col),'*c-','LineWidth',2);hold off;
% % semilogy(Eb_N0_dB,simBer,'mx-','LineWidth',2);
% % 
% ax = gca;
% xlabel('SNR (dB)');
% ax.XAxis.LineWidth = 2;
% ylabel('Bit Error Rate(BER)');
% ax.YAxis.LineWidth = 2;
% grid on;
% legend('Rx1 proposed-sim', 'Rx2 proposed-sim', 'Rx1 proposed', 'Rx2 proposed','EVE-Rx1 proposed', 'EVE-Rx2 proposed', 'Rx1 with Dvsty proposed', 'Rx2 with Dvsty proposed', 'EVE-Rx1 with Dvsty proposed', 'EVE-Rx2 with Dvsty proposed','theory (nTx=1,nRx=1)','sim (nTx=2, nRx=2, ZF)','fontweight','bold');
% %%
% figure(2);
% col=length(d);
% % tx2ber = BER_EVE_Tx2;
% %change line 76 and 99 tohave relationshionship between user channel and
% %EVE
% % semilogy(flip(d(1,1:col)),BER_EVE_Tx1(1,1:col),'*g--','LineWidth',2);hold on
% % semilogy(flip(d(1,1:col)),tx2ber(1,1:col),'ob--','LineWidth',2); hold on;
% % semilogy(flip(d(1,1:col)),BER_EVE_Tx2(1,1:col),'or--','LineWidth',2); hold on;
% % 
% % 
% % ax = gca;
% % xlabel('Distance (m)');
% % ax.XAxis.LineWidth = 2;
% % ylabel('Bit Error Rate (BER)');
% % ax.YAxis.LineWidth = 2;
% % grid on;
% % legend('EVE-Tx1', 'EVE-Tx2','EVE','fontweight','bold')
% 
% %%
% % w=diag(BER_EVE_Tx2);
% % 
% % w(:,1)=BER_EVE_Tx2;
% % w(1,:)=BER_EVE_Tx2;
% % 
% % for i=3:length(BER_EVE_Tx2)
% %     j=2;
% %     while(w(i,j)==0 && j<length(BER_EVE_Tx2))
% %         w(i,j)=w(i,j-1);
% %         j=j+1;
% %     end
% % end
% % 
% % for j=3:length(BER_EVE_Tx2)
% %     i=2;
% %     while(w(i,j)==0 && j<length(BER_EVE_Tx2))
% %         w(i,j)=w(i-1,j);
% %         i=i+1;
% %     end
% % end
% % 
% % bar3(w)
% % %%
% % [x,y]=meshgrid(0:1:80,0:1:80);
% % surf(x,y,w)
% % zlabel('SNR (dB)');
% % ylabel('Bit error rate (BER)');
% % %%
% % plot(abs(x2))
% 
% %% TER
% 
% figure(2);
% semilogy(EbN0,TER_Tx1Rx1,'og-','LineWidth',2);hold on;
% semilogy(EbN0,TER_Tx2Rx2,'*b--','LineWidth',2);hold on;
% % semilogy(EbN0,TER_EVE_Tx1,'vc-','LineWidth',2);hold on;
% % semilogy(EbN0,TER_EVE_Tx2,'^r--','LineWidth',2);hold on; 
% semilogy(EbN0,TER_Tx1Rx1_DIV,'-pm','LineWidth',2);hold on;
% semilogy(EbN0,TER_Tx2Rx2_DIV,'--pk','LineWidth',2);hold on;
% 
% ax = gca;
% xlabel('SNR (dB)');
% ax.XAxis.LineWidth = 2;
% ylabel('Throughput rate');
% ax.YAxis.LineWidth = 2;
% grid on;
% legend('Rx1','Rx2','Rx1 to EVE','Rx2 to EVE','Rx1 with Div','Rx2 with Div','fontweight','bold');
% 
% % %% PER
% % figure(3);
% % semilogy(EbN0(1,1:col),PERu11p(1,1:col),'pg-','LineWidth',2); hold on;
% % semilogy(EbN0(1,1:col),PERu21p(1,1:col),'pb--','LineWidth',2);hold on;
% % semilogy(EbN0(1,1:col),PBERu1E1p(1,1:col),'pc-','LineWidth',2); hold on;
% % semilogy(EbN0(1,1:col),PBERu2E1p(1,1:col),'pr--','LineWidth',2); hold off;
% % 
% % xlabel('SNR (dB)');
% % ax.XAxis.LineWidth = 2;
% % ylabel('Packet Error Rate (PER)');
% % ax.YAxis.LineWidth = 2;
% % grid on;
% % legend('Rx1','Rx2','Rx1 to EVE','Rx2 to EVE','fontweight','bold');
% % 
% % %% PAPR
% % PAPRtx1dB=10*log10(PAPRtx1); % in dB
% % [F1,X1] = ecdf(PAPRtx1dB); %CDF
% % 
% % PAPRtx2dB=10*log10(PAPRtx2); % in dB
% % [F2,X2] = ecdf(PAPRtx2dB); %CDF
% % 
% % PAPRtx3dB=10*log10(PAPRtx3); % in dB
% % [F3,X3] = ecdf(PAPRtx3dB); %CDF
% % 
% % PAPRtx4dB=10*log10(PAPRtx4); % in dB
% % [F4,X4] = ecdf(PAPRtx4dB); %CDF
% % figure(4)
% % semilogy(X1,1-F1,'r','LineWidth',2); hold on; % 1-F1 just to draw the CCDF
% % semilogy(X2,1-F2,'g','LineWidth',2); hold on;
% % semilogy(X3,1-F3,'m','LineWidth',2); hold on; % 1-F1 just to draw the CCDF
% % semilogy(X4,1-F4,'c','LineWidth',2); 
% % xlabel('PAPR');
% % ylabel('CCDF');
% % grid on;
% % legend('Tx1-proposed','Tx2-proposed','Tx1-conventional','Tx2-conventional');
%%
plot(diag(abs(p11)))
hold on
plot(diag(abs(p1)),'g--')
legend('my method','your recomended method')
title('magnitude of p1')


