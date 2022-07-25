warning off;
% clc
clear all

Nbps=1; 
M=2^Nbps;       % Modulation order=2/4/6 for QPSK/16QAM/64QAM
EbN0=1:1:40;    % EbN0
N_iter=100;    % no of packets   % Number of iterations for each EbN0
Nframe=1;       % no. of sybmols in on % Number of symbols per frame
norms=[1 sqrt(2) 0 sqrt(10) 0 sqrt(42)];    % BPSK 4-QAM 16-QAM
NgType=1;       % NgType=1/2 for cyclic prefix/zero padding
if NgType==1, nt='CP';  elseif NgType==2, nt='ZP';   end

% Bob
PowerdB=[0 -8 -17 -21 -25]; % Channel tap power profile 'dB'
Delay=[0 3 5 6 8];          % Channel delay 'sample'
Power=10.^(PowerdB/10);     % Channel tap power profile 'linear scale'
Power=Power/(sum(Power));
Ntap=length(PowerdB);       % Chanel tap number
Lch=Delay(end)+1;           %Channel length

% Eve
PowerdBE=[0 -8 -17 -21 -25];
DelayE=[0 3 5 6 8];
PowerE=10.^(PowerdBE/10);       % Channel tap power profile 'linear scale'
NtapE=length(PowerdBE);         % Chanel tap number
LchE=DelayE(end)+1;             % Channel length(may be diff than no taps)

Nfft=64;                        % FFT size
Ng=64;                          % Ng=0: Guard interval length
Nsym=Nfft+Ng;                   % Symbol duration
Nvc=0;                          % Nvc=0: no virtual carrier
Nused=Nfft; 
%%
bpskMod = comm.BPSKModulator;
bpskDemod = comm.BPSKDemodulator;
%%
 
for i=1:length(EbN0)

 BERu1=0;
 BERu2=0;
 BERu1E=0;
 BER=0;
 
 bittot=0;
 
   for m=1:N_iter    
       
       % Define channels
       %H11
       channel=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);
       h11=zeros(1,Lch); 
       h11(Delay+1)=channel;
       h11r1=h11;
       H11r1=fft([h11r1 zeros(1,Nfft-Lch)]);
     
       %H21
       channel=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);
       h21=zeros(1,Lch); 
       h21(Delay+1)=channel;
       h21r1=h21;
       H21r1=fft([h21r1 zeros(1,Nfft-Lch)]);
      
       %Eve-R1
       channel=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);
       hE1=zeros(1,Lch); 
       hE1(Delay+1)=channel;
       hEr1=hE1;
       HEr1=fft([hEr1 zeros(1,Nfft-Lch)]);
       
       %H12
       channel=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);
       h12=zeros(1,Lch); 
       h12(Delay+1)=channel;              
       h12r2=h12;
       H12r2=fft([h12r2 zeros(1,Nfft-Lch)]);
      
       %H22;
       channel=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);
       h22=zeros(1,Lch); 
       h22(Delay+1)=channel;       
       h22r2=h22;
       H22r2=fft([h22r2 zeros(1,Nfft-Lch)]);
     
       %Eve-R2
       channel=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);
       hE2=zeros(1,Lch); 
       hE2(Delay+1)=channel;
       hEr2=hE2;
       HEr2=fft([hEr2 zeros(1,Nfft-Lch)]);
       
       
       %Data x1 and x2
      X1 = randi([0,M-1], Nused*Nframe,1)';
%       X1mod= bpskMod(X1')/norms(Nbps); %modulate
      X1mod= qammod(X1,M,'gray')/norms(Nbps);
      
      if NgType~=2
          x_GI=zeros(1,Nframe*Nsym);
      elseif NgType==2
          x_GI= zeros(1,Nframe*Nsym+Ng);
      end
     
      X2= randi([0,M-1], Nused*Nframe,1)'; 
%       X2mod= bpskMod(X2')/norms(Nbps);%modulate
      X2mod= qammod(X2,M,'gray')/norms(Nbps);
      if NgType~=2
          x_GI2=zeros(1,Nframe*Nsym);
      elseif NgType==2
          x_GI2= zeros(1,Nframe*Nsym+Ng);
      end


     % Precoder matrices
       A = diag(H11r1);
       B = diag(H12r2);
       C = diag(H21r1);
       D = diag(H22r2);
       
       p1=-D/(C*B-D*A);
       p3=C/(C*B-D*A);

       p2=-B/(D*A-C*B);
       p4=A/(D*A-C*B);
       
  %%    
      %Transmiters
%       Tx1=(1/sqrt(2))*X1mod*(p1.*1/sqrt(abs(p1))) + (1/sqrt(2))*X2mod*(p2.*1/sqrt(abs(p2)));
%       %Antenna 2    
%       Tx2=(1/sqrt(2))*X1mod*(p3.*1/sqrt(abs(p3))) + (1/sqrt(2))*X2mod*(p4.*1/sqrt(abs(p4)));
      p1nom = (norm(diag(p1),2))^2;
      p2nom = (norm(diag(p2),2))^2;
      p3nom = (norm(diag(p3),2))^2;
      p4nom = (norm(diag(p4),2))^2;
      
        vect1 = [diag(p1); diag(p2)];
        vect1nom = norm(vect1,2);
        vect2 = [diag(p3); diag(p4)];
        vect2nom = norm(vect2,2);
      
      Tx1=(1/sqrt(2))*X1mod*(p1*1/vect1nom) + (1/sqrt(2))*X2mod*(p2*1/vect1nom);
      %Antenna 2    
      Tx2=(1/sqrt(2))*X1mod*(p3*1/vect2nom) + (1/sqrt(2))*X2mod*(p4*1/vect2nom);
      
%%
    
      %Transmission
      x1= ifft(Tx1);
%       x1= ifft(Tx1)/norm(ifft(Tx1),2);  %% This is how i am normalizign it :step1
      x_GI1= guard_interval(Ng,Nfft,NgType,x1);       
      y11=conv(x_GI1,h11r1*vect1nom);   %RX1 TX1
      y21=conv(x_GI1,h21r1*vect1nom);   %RX2 TX1
      yE1=conv(x_GI1,hEr1);    %RXE TX1
         
      x2= ifft(Tx2);
%       x2= ifft(Tx2)/norm(ifft(Tx2),2);
      x_GI2= guard_interval(Ng,Nfft,NgType,x2);       
      y12=conv(x_GI2,h12r2*vect2nom);  %RX1 TX2
      y22=conv(x_GI2,h22r2*vect2nom); %RX2 TX2
      yE2=conv(x_GI2,hEr2);  %RXE TX2   

%        snr = ebn0/b, b is spectral efficiency
       snr = EbN0(i)+10*log10(Nbps*(Nused/Nfft));
       
       %Noise 1
%        EEs=sum(x1).^2/length(x1);
       EEs=1;
       noise_mag  = sqrt(0.5*EEs*((10.^(-snr/20))));
       n1 = (1/sqrt(2))*noise_mag*(randn(size(y11))+1i*randn(size(y11)));

       %Noise 2
%        EEs=sum(x2).^2/length(x2);   
       EEs=1;
       noise_mag = sqrt(0.5*EEs*((10.^(-snr/10))));    
       n2 = (1/sqrt(2))*noise_mag*(randn(size(y12))+1i*randn(size(y12)));
       
       %==============Received signals=================================
%        Y1 = y11*norm(ifft(Tx1),2) + y12*norm(ifft(Tx2),2) + n1;
       Y1 = y11 + y12 + n1;       
%        Y2 = y21*norm(ifft(Tx1),2) + y22*norm(ifft(Tx2),2) + n2;
       Y2 = y21 + y22 + n2;
       
       
       Y1new = fft(remove_GI(Ng,Nsym,NgType,Y1));
       Y2new = fft(remove_GI(Ng,Nsym,NgType,Y2));
        
       %Demodulation
%        UU1i=bpskDemod(Y1new'*norms(Nbps));
       UU1i=qamdemod(Y1new*norms(Nbps),M,'gray');
%        UU2i=bpskDemod(Y2new'*norms(Nbps));
       UU2i=qamdemod(Y2new*norms(Nbps),M,'gray');

        
       BERu1= BERu1+sum(sum(de2bi(UU1i,Nbps)~=de2bi(X1,Nbps)));
       BERu2= BERu2+sum(sum(de2bi(UU2i,Nbps)~=de2bi(X2,Nbps)));

       bittot=bittot+Nfft*Nbps; 
   end
 
%BER
 BERu11(i)=BERu1/bittot;
 BERu21(i)=BERu2/bittot;
end
%%
[EbN01,BER_sim]=BERofBPSKinAWGN();
EbN0Lin = 10.^(EbN0/10);
theoryBer = 0.5.*(1-sqrt(EbN0Lin./(EbN0Lin+1)));
% [EbN0,ber]=BERofQPSKinAWGN();
% [EbN0,simBer]=OFDMberINawgn();
%%
figure(1);
semilogy(EbN0,BERu11,'or-','LineWidth',2); hold on;
semilogy(EbN0,BERu21,'*b-','LineWidth',2);hold on;
semilogy(EbN01,BER_sim,'-k*','LineWidth',2);  hold on;
semilogy(EbN0,theoryBer,'-g*','LineWidth',2);  hold on;

% semilogy(EbN0,ber,'--g*','LineWidth',2);  hold on;
% semilogy(EbN0,simBer,'--c*','LineWidth',2);  hold on;
% axis([min(EbN0) max(EbN0) 10^(-5) 1]);
xlabel('SNR (dB)');
ylabel('Bit Error Rate(BER)');
legend('TX1','TX2', 'Rayleigh-Simulation','fontweight','bold');

% legend('TX1','TX2', 'BPSK AWGN perfom','QPSK in AWGN perfom','BPSKofdmAWGN','MIMOawgn','fontweight','bold');
grid on;
%%
% x12=(diag(H11r1)*p2)+(diag(H12r2)*p4);
% 
% newNoise= diag(x21);
% newNoise=newNoise';
% noise=[];
% for i=1:12500
%     noise=[noise newNoise];
% end
% 

