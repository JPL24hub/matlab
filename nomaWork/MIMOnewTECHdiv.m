
% function [nErr1, nErr2] = MIMOnewTECHdiv(iter)
clear all

Eb_N0_dB = [-3:35]; % multiple Eb/N0 values
n_iter = 2000;
symbol_len=64;

nTx = 2;
nRx = 2;

for ii = 1:length(Eb_N0_dB)
    BER1=0;
    BER2=0;
    bittot=0;
    
    for k = 1:n_iter
        
        % channels
        h11 = 1/sqrt(2)*[randn(symbol_len,1) + 1i*randn(symbol_len,1)]; % Rayleigh channel
        h12 = 1/sqrt(2)*[randn(symbol_len,1) + 1i*randn(symbol_len,1)]; % Rayleigh channel
        h21 = 1/sqrt(2)*[randn(symbol_len,1) + 1i*randn(symbol_len,1)]; % Rayleigh channel
        h22 = 1/sqrt(2)*[randn(symbol_len,1) + 1i*randn(symbol_len,1)]; % Rayleigh channel
        
        x1 = rand(symbol_len,1)>0.5; % generating 0,1 with equal probability
        s1 = 2*x1-1; % BPSK modulation 0 -> -1; 1 -> 0
        
        x2 = rand(symbol_len,1)>0.5; % generating 0,1 with equal probability
        s2 = 2*x2-1; % BPSK modulation 0 -> -1; 1 -> 0
        
         % Precoder matrices
        A = diag(h11);
        B = diag(h12);
        C = diag(h21);
        D = diag(h22);
       
        AA = abs(A)^2;
        BB = abs(B)^2;
        CC = abs(C)^2;
        DD = abs(D)^2;
       
       
        p1=(-D*(AA+BB))/(C*B-D*A);
        p3=(C*(AA+BB))/(C*B-D*A);

        p2=(-B*(CC+DD))/(D*A-C*B);
        p4=(A*(CC+DD))/(D*A-C*B); 
        
        
        vect1 = [diag(p1); diag(p2)];
        vect1nom = norm(vect1,2);
        vect2 = [diag(p3); diag(p4)];
        vect2nom = norm(vect2,2);

%         
%         u1=(1/sqrt(2))*diag(s1)*(p1*1/vect1nom) + (1/sqrt(2))*diag(s2)*(p2*1/vect1nom);
%         %Antenna 2    
%         u2=(1/sqrt(2))*diag(s1)*(p3*1/(vect2nom)) + (1/sqrt(2))*diag(s2)*(p4*1/(vect2nom));
        
        % Transmitters
        u1 = diag(s1)*(p1) + diag(s2)*(p2);   
        u2 = diag(s1)*(p3) + diag(s2)*(p4);
        
        % noise
        n1 = 1/sqrt(2)*[randn(symbol_len,1) + 1i*randn(symbol_len,1)];  % Rayleigh channel
        n1 = n1* 10^(-Eb_N0_dB(ii)/20);
        
        n2 = 1/sqrt(2)*[randn(symbol_len,1) + 1i*randn(symbol_len,1)];  % Rayleigh channel
        n2 = n2* 10^(-Eb_N0_dB(ii)/20);
        
        % reception
%         y1 = diag(h11)*diag(u1)*vect1nom + diag(h12)*diag(u2)*vect2nom + diag(n1);
%         y2 = diag(h21)*diag(u1)*vect1nom + diag(h22)*diag(u2)*vect2nom + diag(n2);
        
                % reception
        y1 = diag(h11)*diag(u1) + diag(h12)*diag(u2) + diag(n1);
        y2 = diag(h21)*diag(u1) + diag(h22)*diag(u2) + diag(n2);
        
        % receiver - hard decision decoding
        x1Hat = real(diag(y1))>0;
        x2Hat = real(diag(y2))>0;
        
        % error calculation
        BER1 = BER1 + size(find([x1 - x1Hat]),1);
        BER2 = BER2 + size(find([x2 - x2Hat]),1);
        
        bittot = bittot + symbol_len;
        
    end
    nErr1(ii) = BER1/bittot;
    nErr2(ii) = BER2/bittot;
end
%%
EbN0Lin = 10.^(Eb_N0_dB/10);
theoryBer = 0.5.*(1-sqrt(EbN0Lin./(EbN0Lin+1)));
semilogy(Eb_N0_dB,nErr1,'om-','LineWidth',2); hold on
semilogy(Eb_N0_dB,nErr2,'*g--','LineWidth',2); hold on
semilogy(Eb_N0_dB,theoryBer,'cd--','LineWidth',2);
grid on
legend('TX1-ideal', 'TX2-ideal', 'Rayleigh-Simulation')
% end

% http://www.dsplog.com/2009/04/13/transmit-beamforming/
% http://www.dsplog.com/2008/08/19/receive-diversity-in-awgn/
% http://www.dsplog.com/2009/03/15/alamouti-stbc-2-receive-antenna/








