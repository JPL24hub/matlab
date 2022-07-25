
  
clear all;
clc
warning off;

Nbps=1; 
M=2^Nbps;           % Modulation order=2/4/6 for QPSK/16QAM/64QAM
EbN0=0:5:40;      % EbN0

N_iter=100;         %1e5;   % no of packets   % Number of iterations for each EbN0
Nframe=1;           % no. of sybmols in on % Number of symbols per frame

norms=[1 sqrt(2) 0 sqrt(10) 0 sqrt(42)];     % BPSK 4-QAM 16-QAM
NgType=1;           % NgType=1/2 for cyclic prefix/zero padding


%% Bob
PowerdB=[0 -8 -17 -21 -25];         % Channel tap power profile 'dB'
Delay=[0 3 5 6 8];                  % Channel delay 'sample'
Power=10.^(PowerdB/10);             % Channel tap power profile 'linear scale'
Power=Power/(sum(Power));
Ntap=length(PowerdB);               % Chanel tap number
Lch=Delay(end)+1;                   %Channel length

%% Eve
PowerdBE=[0 -8 -17 -21 -25];
DelayE=[0 3 5 6 8];
PowerE=10.^(PowerdBE/10);           % Channel tap power profile 'linear scale'
NtapE=length(PowerdBE);             % Chanel tap number
LchE=DelayE(end)+1;                 % Channel length(may be diff than no taps)
%%

Nfft=64;                            % FFT size
Ng=16;                              % Ng=3; Ng=0: Guard interval length

Nsym=Nfft+Ng;                       % Symbol duration
Nvc=0;                              % Nfft/4; Nvc=0: no virtual carrier
Nused=Nfft ;                        % Nvc;
 
for i=1:length(EbN0)

                                    % Neb=0; Ntb=0; % Initialize the number of error/total bits
 BERu1=0;
 BERu2=0;
 BERu1E=0;
 BERu2E=0;
 
 bittot=0;
 
   for m=1:N_iter
     %% Round 1 
      %user1-atenna1
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
      

       %% Round 2 -Antenna 2
        %h_ab  a=usr b=antena
      %user1  
       channel=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);%%%%%
       h12=zeros(1,Lch); 
       h12(Delay+1)=channel; % cir: channel impulse response
       h12r2=h12;
       H12r2=fft([h12r2 zeros(1,Nfft-Lch)]); % Channel frequency response
      %user2
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
      

       %%
      % User1 Tx______________________________________________________________
      X= randi([0,M-1], Nused*Nframe,1)'; % bit: integer vector
      Xmod= qammod(X,M,'gray')/norms(Nbps); 
      
      if NgType~=2, x_GI=zeros(1,Nframe*Nsym);
      elseif NgType==2, x_GI= zeros(1,Nframe*Nsym+Ng);
        % Extend an OFDM symbol by Ng zeros 
      end
      
      % User2 Tx______________________________________________________________
    
      X2= randi([0,M-1], Nused*Nframe,1)'; % bit: integer vector
      Xmod2= qammod(X2,M,'gray')/norms(Nbps);
      
      if NgType~=2, x_GI2=zeros(1,Nframe*Nsym);
      elseif NgType==2, x_GI2= zeros(1,Nframe*Nsym+Ng);
        % Extend an OFDM symbol by Ng zeros 
      end

        %Modulted signal for user1 and user2      
      U1=Xmod;
      U2=Xmod2;
        %% calculate r1 r2
 
        A = diag(H11r1);
        B = diag(H12r2);
        C = diag(H21r1);
        D = diag(H22r2);
        
        S = ((A+B)*diag(U2));
        T = ((C+D)*diag(U1));
        
%         S = diag(U1);
%         T = diag(U2);
        
        r2 = ((S*C)-(A*T))/((D*A)-(C*B));
        r1 = (-S-(B*r2))/A;
        
        %WP2 MRC
%         A = abs(diag(H11r1)).^(2);
%         B = abs(diag(H12r2)).^(2);
%         C = abs(diag(H21r1)).^(2);
%         D = abs(diag(H22r2)).^(2);
%         
%         S = -((A+B)*diag(U2));
%         T = -((C+D)*diag(U1));
%         
%         r2 = ((S*C)-(A*T))/((B*C)-(A*D));
%         r1 = (T-(D*r2))/C;
        

        %% superposition
       Tx1=sqrt(1)*U1+sqrt(1)*U2+sqrt(1)*diag(r1)';
       Tx2=sqrt(1)*U1+sqrt(1)*U2+sqrt(1)*diag(r2)';
        
        %% for user 1. normalize and add guard intervals
        
         x1= ifft(Tx1)/norm(ifft(Tx1),2);  %% This is how i am normalizign it :step1
         norm(ifft(Tx1),2);

         x_GI1= guard_interval(Ng,Nfft,NgType,x1);       
         y11= conv(x_GI1,h11r1); %u1 r1
         y21=(conv(x_GI1,h21r1)); %u2 r1
         yE1= conv(x_GI1,hEr1); %u2 r1

         snr = EbN0(i)+10*log10(Nbps*(Nused/Nfft)); % SNR vs. Eb/N0
         EEs=sum(x1).^2/length(x1);    
         noise_mag  = sqrt(0.5*EEs*((10.^(-snr/10))));
         N1 = noise_mag*(randn(size(y11))+1i*randn(size(y11)));
         
         y_G11r1 = norm(ifft(Tx1),2)*y11+ N1;  %% This is how i am normalizign : step2
         y_G21r1 = norm(ifft(Tx1),2)*y21+ noise_mag*(randn(size(y21))+1i*randn(size(y21)));
         y_G21E1 = norm(ifft(Tx1),2)*yE1+ noise_mag*(randn(size(yE1))+1i*randn(size(yE1)));
 
         
         Y11= fft(remove_GI(Ng,Nsym,NgType,y_G11r1));
         Y21= fft(remove_GI(Ng,Nsym,NgType,y_G21r1));
         YE1= fft(remove_GI(Ng,Nsym,NgType,y_G21E1));
 
         %% for user 2. normalize and add guard intervals
         
        x2= ifft(Tx2)/norm(ifft(Tx2),2);
        
        x_GI2= guard_interval(Ng,Nfft,NgType,x2);       
        y12= conv(x_GI2,h12r2); %u1 r2
        y22= conv(x_GI2,h22r2); %u2 r2
        yE2= conv(x_GI2,hEr2); %u2 r1    
  
        EEs=sum(x2).^2/length(x2);    
        noise_mag  = sqrt(0.5*EEs*((10.^(-snr/10))));
        N2=noise_mag*(randn(size(y12))+1i*randn(size(y12)));
        
        y_G12r2 = norm(ifft(Tx2),2)* y12+ N2;
        y_G22r2 = norm(ifft(Tx2),2)* y22+ N2; %noise_mag*(randn(size(y22))+1i*randn(size(y22)));
        y_G21E2 = norm(ifft(Tx2),2)*yE2+ N2; %noise_mag*(randn(size(yE2))+1i*randn(size(yE2)));
 
        Y12= fft(remove_GI(Ng,Nsym,NgType,y_G12r2));
        Y22= fft(remove_GI(Ng,Nsym,NgType,y_G22r2));
        YE2= fft(remove_GI(Ng,Nsym,NgType,y_G21E2));
         
         %%
%          YU1= (H11r1)'.*transpose(Y11)+(H12r2)'.*transpose(Y12);
%          YU2= (H21r1)'.*transpose(Y21)+(H22r2)'.*transpose(Y22);
%          YUe= YE1+YE2;
         
         YU1= Y11+Y12;
         YU1 = YU1/(A+B);
         YU2= Y21+Y22;
         YU2 =YU2/(C+D);
         YUe= YE1+YE2;


        UU1i=qamdemod(YU1*norms(Nbps),M,'gray');
        UU2i=qamdemod(YU2*norms(Nbps),M,'gray');
        UUee=qamdemod(YUe*norms(Nbps),M,'gray');
  

        %%
        BERu1=BERu1+sum(sum(de2bi(UU1i,Nbps)~=de2bi(X,Nbps)));
        BERu2= BERu2+sum(sum(de2bi(UU2i,Nbps)~=de2bi(X2,Nbps)));

        %% External Eve
        BERu1E=BERu1E+sum(sum(de2bi(UUee,Nbps)~=de2bi(X,Nbps)));
        BERu2E=BERu2E+sum(sum(de2bi(UUee,Nbps)~=de2bi(X2,Nbps)));
 

        bittot=bittot+Nfft*Nbps; 

   end 
 
%BER
 BERu11(i)=BERu1/bittot;
 BERu21(i)=BERu2/bittot;
 BERu1E1(i)=BERu1E/bittot;
 BERu2E1(i)=BERu2E/bittot;

   
end
% 
figure(1);
semilogy(EbN0,BERu11,'or-','LineWidth',2); hold on;
semilogy(EbN0,BERu21,'*b-','LineWidth',2);hold on;
semilogy(EbN0,BERu1E1,'vc-','LineWidth',2); hold on;
semilogy(EbN0,BERu2E1,'^r-','LineWidth',2); 
xlabel('SNR (dB)');
ylabel('Bit Error Rate(BER)');
legend('BER(Bob1)','BER(Bob2)','Eve listento Bob1','Eve listento Bob2');
