%% Analysis Secrecy Capacity of dual hop DF Relaying system %%
% Elham Nosrati, March 2014
% S(-500,0), R(0,0), D(500,0), Eve(0,1000) meters
% Path loss coefficient: alfa=3.5
% I have used _db for the SNR in dB
% for simple ransmission Standard deviation of each channel is 1/d
% SD_rd=1/d_rd , ch = sqrt(variance)*( randn(1,N)+j*randn(1,N));
%% Location of nodes and Power allocation %%%
clear all
clc
N= 10^3; % Monte Carlo
Pt_dBm=[-5:1:35];  %Global Transit power in dBm  
for l=1:length(Pt_dBm)
Ptotal(l)=(10^(Pt_dBm(l)/10))/1000;  % Global Transmit power in Watt
end
    
Rs=0.1;   % Rs=Secrecy rate
for i=1:length(Ptotal)
    Pt=Ptotal(i);% Total tranmited Power
    Ps=Pt/2;% Source transmit power
    Pr=Pt/2; % Relay transmit power
    
    alfa=3.5; % Path loss exponent
    % Source, Relay and Destination locations
    d_sr=500;  %distance between the source and the relay
    d_se=1000*sqrt(1.25);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Legitimate Channel%%%%%%%%
    %% Fading of legitimate channel
        channel_SR=sqrt(1/d_sr).*[randn(1,N)+j*randn(1,N)];% Channel coefficients for S->R Link           
             
        Amp_CH_SR=abs(channel_SR); % channel amplitude=|h_SR|
        CH_SR=Amp_CH_SR.^2; % |h_SR|^2
        
       
 %% AWGN in legitimate channel    
%         n_SR=(1/sqrt(2)).*(randn(1,N)+j*randn(1,N));%AWGN at Relay for S->R Link
%         n_RD=(1/sqrt(2)).*(randn(1,N)+j*randn(1,N));%AWGN at Desination  for R->D Link
%         
%         N_SR=var(n_SR);   % Noise Variance S->R Link
%         N_RD=var(n_RD);   % Noise Variance R->D Link
 
% we consider noise variance is (-60) dBm at all nodes
        No=(10^(-60/10))/1000;
        N_SR=No;   % Noise Variance S->R Link
  
%% Received SNR at Legitimate nodes
    
                %First hop SNR at Relay
                SNR1_linear=(Ps)*(1/(d_sr^alfa))*(CH_SR/N_SR);
                snr1_db=10*log10(SNR1_linear);
                
                SNR_SR=mean(SNR1_linear);
                SNR_SR_db=10*log10(SNR_SR);
                
               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Wiretap channel%%%%%%%%%%%%%%%%%%%%     
%% Fading of Wiretap channel
        channel_SE=(1/ d_se).*[randn(1,N)+j*randn(1,N)];% Channel coefficients for S->E Link            
              
        Amp_CH_SE=abs(channel_SE); %   channel amplitude=|h_SR|
        CH_SE=Amp_CH_SE.^2; % |h_SR|^2
       
 %% AWGN in wiretap channel  
%         n_SE=(1/sqrt(2)).*(randn(1,N)+j*randn(1,N));%AWGN at Relay for S->E Link
%         n_RE=(1/sqrt(2)).*(randn(1,N)+j*randn(1,N));%AWGN at Desination  for R->E Link
%         
%         N_SE=var(n_SE);   % Noise Variance S->R Link
%         N_RE=var(n_RE);   % Noise Variance R->D Link
        
% we consider noise variance is (-60) dBm at all nodes
        N_SE=No;   % Noise Variance S->R Link
     
        
 %% Received SNR at Eavesdropper
                %First hop SNR at Eve
                SNR_E1_linear=(Ps)*(1/(d_se^alfa))*(CH_SE/N_SE);
                snr_E1_db=10*log10(SNR_E1_linear);
                                
                SNR_SE=mean(SNR_E1_linear);
                SNR_SE_db=10*log10(SNR_SE);
                
                
               % Total Received at Eve
               SNR_E_linear=SNR_SE                     
               snr_E_db=10*log10(SNR_E_linear); 
               
   
%% Calculating secrecy capacity 
      
% Secrecy Capacity of the system
% Cs Analysis
   
Cs_Analys(i)= (1/2)*log2((1+SNR_SR)/(1+SNR_E_linear));
end
%% Ploting Results
% Ploting Cs
plot(Pt_dBm,Cs_Analys,'-ro');
xlabel('Global Transit power[dBm]');
ylabel('Secrecy Capacity');
    
    
    
