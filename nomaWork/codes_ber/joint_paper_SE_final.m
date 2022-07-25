%% ----------- necessary definitions, channels etc -------------------
clear all;
warning off;

%% parameters
%Important parameters: L , Modulation, chanel coding , fft size, frmae size 
Nbps=1; 
M=2^Nbps;  % Modulation order=2/4/6 for QPSK/16QAM/64QAM
EbN0 = [0:5:30];    % EbN0
N_iter=100;   % no of packets   % Number of iterations for each EbN0
Nframe=1;         % no. of sybmols in on % Number of symbols per frame
sigPow=0;         % Signal power initialization
norms=[1 sqrt(2) 0 sqrt(10) 0 sqrt(42)];     % BPSK 4-QAM 16-QAM

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

%% ---- NOMA new
% N = 10^5;
N = 64;
d1 = 1000; d2 = 500;    %Distances of users from base station (BS)
a1 = 0.75; a2 = 0.25;   %Power allocation factors
eta = 4;                %Path loss exponent

%Generate rayleigh fading coefficient for both users
h1 = sqrt(d1^-eta)*(randn(1,N)+1i*randn(1,N))/sqrt(2);
h2 = sqrt(d2^-eta)*(randn(1,N)+1i*randn(1,N))/sqrt(2);

g1 = (abs(h1)).^2;
g2 = (abs(h2)).^2;

Pt = 0:5:30;                %Transmit power in dBm
pt = (10^-3)*10.^(Pt/10);   %Transmit power in linear scale
BW = 10^6;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)
p = length(Pt);

for i=1:length(EbN0)
%% ------------ user channels --------------------------------
%% user 1
       channel=(randn(1,Ntap)+j*randn(1,Ntap)).*sqrt(Power/2);%%%%%
       h11=zeros(1,Lch); 
       h11(Delay+1)=channel; % cir: channel impulse response
       h11A1=h11;
       H11A1=fft([h11A1 zeros(1,Nfft-Lch)]); % Channel frequency response
       
       channel=(randn(1,Ntap)+j*randn(1,Ntap)).*sqrt(Power/2);%%%%%
       h12=zeros(1,Lch); 
       h12(Delay+1)=channel; % cir: channel impulse response
       h12A2=h12;
       H12A2=fft([h12A2 zeros(1,Nfft-Lch)]); % Channel frequency response
%% user 2
       channel=(randn(1,Ntap)+j*randn(1,Ntap)).*sqrt(Power/2);%%%%%
       h21=zeros(1,Lch); 
       h21(Delay+1)=channel; % cir: channel impulse response
       h21A1=h21;
       H21A1=fft([h21A1 zeros(1,Nfft-Lch)]); % Channel frequency response

       channel=(randn(1,Ntap)+j*randn(1,Ntap)).*sqrt(Power/2);%%%%%
       h22=zeros(1,Lch); 
       h22(Delay+1)=channel; % cir: channel impulse response
       h22A2=h22;
       H22A2=fft([h22A2 zeros(1,Nfft-Lch)]); % Channel frequency response
       
       %% normalize the channels to simulate the doewlink NOMA in aldababsa
       H11A1  = H11A1./norm(H11A1);
       H21A1 = H21A1.*(10/norm(H21A1));
       H12A2 = H12A2./(norm(H12A2));
       H22A2 = H22A2./(10/norm(H22A2));
%% modulated symbols
X1= randi([0,M-1], Nused*Nframe,1)'; % 0 1 2 3
Xmod1= qammod(X1,M,'gray')/norms(Nbps); 
X2= randi([0,M-1], Nused*Nframe,1)'; 
Xmod2= qammod(X2,M,'gray')/norms(Nbps);

       %% ---------- noise generation --------------
% i = 1;
snr = EbN0(i)+10*log10(Nbps*(Nused/Nfft)); 
EEs=sum(X1).^2/length(X1);    
noise_mag  = sqrt(0.5*EEs*((10.^(-snr/10))));
user1_n1 = noise_mag*(randn(size(X1))+j*randn(size(X1)));
EEs2=sum(X2).^2/length(X2);    
noise_mag2  = sqrt(0.5*EEs2*((10.^(-snr/10))));
user2_n2 = noise_mag2*(randn(size(X2))+j*randn(size(X2)));
% n1 = user1_n1;
% n2 = user2_n2;
%% --------------------------- WP1 ---------------------------------------------
 M_1a_D1=abs(diag(H22A2)).^(2)*(abs(diag(H22A2)).^(2)*abs(diag(H11A1)).^(2)-abs(diag(H21A1)).^(2)*abs(diag(H12A2)).^(2))^-1;
 M_1b_D1=-(M_1a_D1*abs(diag(H21A1)).^(2))/abs(diag(H22A2).^(2));

 M_2a_D1=abs(diag(H12A2)).^(2)*(abs(diag(H21A1)).^(2)*abs(diag(H12A2)).^(2)-abs(diag(H22A2)).^(2)*abs(diag(H11A1)).^(2))^-1;
 M_2b_D1=-(M_2a_D1*abs(diag(H11A1)).^(2))/abs(diag(H12A2)).^(2);

 %  y1Hat_D1 = (abs(diag(H11A1)).^(2)*diag(M_1a_D1) + abs(diag(H12A2)).^(2)*diag(M_1b_D1))*Xmod1 + (abs(diag(H11A1)).^(2)*diag(M_2a_D1) + abs(diag(H12A2)).^(2)*diag(M_2b_D1))*Xmod2 +n1;
%  y2Hat_D1 = (abs(diag(H21A1)).^(2)*diag(M_1a_D1) + abs(diag(H22A2)).^(2)*diag(M_1b_D1))*Xmod1 + (abs(diag(H21A1)).^(2)*diag(M_2a_D1) + abs(diag(H22A2)).^(2)*diag(M_2b_D1))*Xmod2 + n2;
  S13 = norm((abs(diag(H11A1)).^(2)*diag(M_1a_D1) + abs(diag(H12A2)).^(2)*diag(M_1b_D1)).*(Xmod1)).^2;
  I13 = (norm((abs(diag(H21A1)).^(2)*diag(M_1a_D1) + abs(diag(H22A2)).^(2)*diag(M_1b_D1)).*Xmod1)).^2;
  n13 = user1_n1;
 
  S23 = (norm((abs(diag(H21A1)).^(2)*diag(M_2a_D1) + abs(diag(H22A2)).^(2)*diag(M_2b_D1)).*Xmod2)).^2;
  I23 =  (norm((abs(diag(H11A1)).^(2)*diag(M_2a_D1) + abs(diag(H12A2)).^(2)*diag(M_2b_D1)).*Xmod2)).^2;
  n23 = user2_n2;
 
  SINR13(i) = S13/(I13 + var(n13));
  SINR23(i) = S23/(I23 + var(n23));
 
  RD1(i) = log2((1+SINR13(i))*(1+SINR23(i)));
  
  %% --------------------------- WP2 (in the paper) ---------------------------------------------

  %% precoding values
 M_1a_D3=diag(H22A2)*(diag(H22A2)*diag(H11A1)-diag(H21A1)*diag(H12A2))^-1;
 M_1b_D3=-(M_1a_D3*diag(H21A1))/diag(H22A2);

 M_2a_D3=diag(H12A2)*(diag(H21A1)*diag(H12A2)-diag(H22A2)*diag(H11A1))^-1;
 M_2b_D3=-(M_2a_D3*diag(H11A1))/diag(H12A2);

%% received signals
%  y1Hat_D3 = (H11A1*diag(M_1a_D3) + H12A2*diag(M_1b_D3))*Xmod1 + (H11A1*diag(M_2a_D3) + H12A2*diag(M_2b_D3))*Xmod2 +n1;
%  y2Hat_D3 = (H21A1*diag(M_1a_D3) + H22A2*diag(M_1b_D3))*Xmod1 + (H21A1*diag(M_2a_D3) + H22A2*diag(M_2b_D3))*Xmod2 + n2;
%  
  S13_D3 = norm(H11A1*diag(M_1a_D3) + H12A2*diag(M_1b_D3).*Xmod1).^2;
 I13_D3 = norm(H21A1*diag(M_1a_D3) + H22A2*diag(M_1b_D3).*Xmod1).^2;
 n13_D3 = user1_n1;

  S23_D3 = norm(H21A1*diag(M_2a_D3) + H22A2*diag(M_2b_D3).*Xmod2).^2;
 I23_D3 = norm(H11A1*diag(M_2a_D3) + H12A2*diag(M_2b_D3).*Xmod2).^2;
 n23_D3 = user2_n2;
 SINR13_D3(i) = S13_D3/(I13_D3 + var(n13_D3));
 SINR23_D3(i) = S23_D3/(I23_D3 + var(n23_D3));
 
 RD3(i) = log2((1+SINR13_D3(i))*(1+SINR23_D3(i)));
%% --------------------------- WP3 ---------------------------------------------
        A = abs(diag(H11A1)).^(2);
        B = abs(diag(H12A2)).^(2);
        C = abs(diag(H21A1)).^(2);
        D = abs(diag(H22A2)).^(2);
        
        S = -((A+B)*diag(Xmod2));
        T = -((C+D)*diag(Xmod1));
        
        r2 = ((S*C)-(A*T))/((B*C)-(A*D));
        r1 = (T-(D*r2))/C;
  
        
   S12 = norm((abs(H11A1).^2 + abs(H12A2).^2).*Xmod1).^2;
   I12 = norm((abs(H21A1).^2 + abs(H22A2).^2).*Xmod1).^2;
   n12 = user1_n1;

   S22 = norm((abs(H21A1).^2 + abs(H22A2).^2).*Xmod2).^2;
   I22 = norm((abs(H11A1).^2 + abs(H12A2).^2).*Xmod2).^2;
   n22 = user2_n2;
 SINR12(i) = S12/(I12 + var(n12));
 SINR22(i) = S22/(I22 + var(n22));
 
 RD2(i) = log2((1+SINR12(i))*(1+SINR22(i)));
 %% --------------------------- WP4 ---------------------------------------------

 d_H11A1 = diag(H11A1); d_H21A1 = diag(H21A1);
d_H12A2 = diag(H12A2); d_H22A2 = diag(H22A2);
bigX1 = -(d_H21A1 + d_H22A2)*diag(Xmod1); bigX2 = -(d_H11A1 + d_H12A2)*diag(Xmod2);
r2 = ((d_H22A2 - d_H21A1*(d_H11A1^-1)*d_H12A2)^-1)*(bigX1 - d_H21A1*(d_H11A1^-1)*bigX2);
r1 = (d_H11A1^-1)*(bigX2 - d_H12A2*r2);

  S14 = norm((H11A1 + H12A2).*Xmod1).^2;
  I14 = norm((H21A1 + H22A2).*(Xmod1)+ (H11A1.*r1 + H12A2.*r2)).^2;
  n14 = user1_n1;

  S24 = norm((H21A1 + H22A2).*Xmod2).^2;
  I24 = norm((H11A1 + H12A2).*(Xmod2)+ (H11A1.*r1 + H12A2.*r2)).^2;
  n24 = user2_n2;
 SINR14(i) = S14/(I14 + var(n14));
 SINR24(i) = S24/(I24 + var(n24));
 
 RD4(i) = log2((1+SINR14(i))*(1+SINR24(i)));
%% --------------- downlink NOMA rate --------------------------------------------------
% H11A1_NOMA  = H11A1./norm(H11A1);
% H21A1_NOMA = H21A1.*(10/norm(H21A1));
% SINR1N(i) = 0.6*(norm(H11A1_NOMA)).^2/(0.4*(norm(H11A1_NOMA)).^2 + var(user1_n1)); 
% SINR2N(i) = 0.4*(norm(H21A1_NOMA)).^2/( var(user2_n2));
% 
% RNOMA(i) = log2((1+SINR1N(i))*(1+SINR2N(i)));
end
%% ---- integrate NOMA from the other code -------
for u = 1:p
    %Calculate SNRs
    gamma_1 = a1*pt(u)*g1./(a2*pt(u)*g1+no);
%     gamma_12 = a1*pt(u)*g2./(a2*pt(u)*g2+no);
    gamma_2 = a2*pt(u)*g2/no;
    
    %Calculate achievable rates
    R1 = log2(1+gamma_1);
%     R12 = log2(1+gamma_12);
    R2 = log2(1+gamma_2);
    
    %Find average of achievable rates
    R1_av(u) = mean(R1);
%     R12_av(u) = mean(R12);
    R2_av(u) = mean(R2);
    
    R_av(u) = R1_av(u) + R2_av(u);
end



plot(EbN0,RD1,'LineWidth',2);
hold on;
plot(EbN0,RD3,'LineWidth',2);
hold on;
plot(EbN0,RD2,'LineWidth',2);
hold on;
plot(EbN0,RD4,'LineWidth',2); 
hold on;
plot(EbN0,R_av,'LineWidth',2); 
xlabel('SNR (dB)');
ylabel('achievable Rate');
% ylim([-10 20]);
legend('CT_{1}','CT_{2}','CT_{3}','CT_{4}','Downlink power-domain NOMA' );
grid on;