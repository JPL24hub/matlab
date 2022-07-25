
clc;clear variables; close all;

N = 10^4; %Number of data samples

EbN0dB = 0:2:40; %Eb/N0 in dB overwhich the performance has to be simulated 

totPower=1; %Total power of LOS path & scattered paths 

K=[1 10]; %A list of Ricial K factors to simulate 

%data generation 
d1=rand(1,N)>0.5; 
d2=rand(1,N)>0.5; 

%BPSK modulation 
xp=2*d1-1; 
xq=2*d2-1; 

%Power allocation factors
alpha=0.2; beta=0.8; 
Ps=1; 

%Do superpostion coding
x=sqrt(alpha*Ps)*xp+sqrt(beta*Ps)*xq; 

%Noise samples (AWGN) with zero mean and unit variance
noise = (randn(1,N)+1i*randn(1,N))/sqrt(2);

%Path loss calculation
Lp = ((3*10^8)/(4*3.14*2*10^9*100))^2;  %Free space loss
Gs = 52.1;                             
beamgain = 100;                     
m=Lp*db2pow(Gs)*beamgain;               %Total loss

%Buffers to hold BERs of user1 and user 2 (One row for every value of K)
simBER_ricean1 = zeros(length(K),length(EbN0dB)); 
simBER_ricean2 = zeros(length(K),length(EbN0dB));

for u = 1:length(K)
    
    %Ricean parameters
    k = K(u);
    s = sqrt(k/(k+1)*totPower); 
    sigma=totPower/sqrt(2*(k+1)); 
    
    %Generate Ricean samples
    h = ((sigma*randn(1,N)+s)+1i*(randn(1,N)*sigma+0)); 


    for i=1:length(EbN0dB) 

        n = noise*10^(-EbN0dB(i)/20); %Scale noise according to required Eb/No
        
        y_ricean=sqrt(m).*h.*x+n;             %Received signal
        y_ricean_cap = y_ricean./(sqrt(m)*h); %Equalize
        
        %User 2 (Direct decoding) (Since beta > alpha)
        r_ricean2=real(y_ricean_cap)>0;         %Direct BPSK demodulation of user 1's data

        %User 1;
        r_ricean_remod = 2*r_ricean2-1;         %Remodulate user 2's data in order to do SIC
        rem = y_ricean_cap - sqrt(beta)*r_ricean_remod; %perform SIC
        r_ricean1 = real(rem)>0;                %Demodulate to get user 1's data
        
        %BER calculation
        simBER_ricean1(u,i) = biterr(r_ricean1,d1)/N;   %BER of user 1
        simBER_ricean2(u,i) = biterr(r_ricean2,d2)/N;   %BER of user 2
    end
    
    %Plot BERs
    semilogy(EbN0dB,simBER_ricean1(u,:),'linewidth',2); hold on; grid on;   
    semilogy(EbN0dB,simBER_ricean2(u,:),'--','linewidth',2);
end
xlabel('E_b/N_o (dB)');
ylabel('BER');
legend('User 1, K = 1','User 2, K = 1','User 1, K = 10','User 2, K = 10');




