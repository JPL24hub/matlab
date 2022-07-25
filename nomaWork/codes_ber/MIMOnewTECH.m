
% clear all

Eb_N0_dB = [-3:35]; % multiple Eb/N0 values
n_iter = 4000;
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
       
        p1=-D/(C*B-D*A);
        p3=C/(C*B-D*A);

        p2=-B/(D*A-C*B);
        p4=A/(D*A-C*B);  
        
        
        vect1 = [diag(p1); diag(p2)];
        vect1nom = norm(vect1,2);
        vect2 = [diag(p3); diag(p4)];
        vect2nom = norm(vect2,2);

        
        u1=(1/sqrt(2))*diag(s1)*(p1*1/vect1nom) + (1/sqrt(2))*diag(s2)*(p2*1/vect1nom);
        %Antenna 2    
        u2=(1/sqrt(2))*diag(s1)*(p3*1/(vect2nom)) + (1/sqrt(2))*diag(s2)*(p4*1/(vect2nom));
        
%         % Transmitters
%         u1 = diag(s1)*(p1) + diag(s2)*(p2);   
%         u2 = diag(s1)*(p3) + diag(s2)*(p4);
        
        % noise
        n1 = 1/sqrt(2)*[randn(symbol_len,1) + 1i*randn(symbol_len,1)];  % Rayleigh channel
        n1 = n1* 10^(-Eb_N0_dB(ii)/20);
        
        n2 = 1/sqrt(2)*[randn(symbol_len,1) + 1i*randn(symbol_len,1)];  % Rayleigh channel
        n2 = n2* 10^(-Eb_N0_dB(ii)/20);
        
        % reception
        y1 = diag(h11)*diag(u1)*vect1nom + diag(h12)*diag(u2)*vect2nom + diag(n1);
        y2 = diag(h21)*diag(u1)*vect1nom + diag(h22)*diag(u2)*vect2nom + diag(n2);
        
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

EbN0Lin = 10.^(Eb_N0_dB/10);
theoryBer = 0.5.*(1-sqrt(EbN0Lin./(EbN0Lin+1)));
theoryBerAWGN = 0.5*erfc(sqrt(10.^(Eb_N0_dB/10))); % theoretical ber
[Eb_N0_dBzf,theoryBer_nRx1,theoryBerMRC_nRx2,simBer]=MIMOwithZF();
[BER_Tx1Rx1_DIV, BER_Tx2Rx2_DIV] = MIMOnewTECHdiv(n_iter);
%%

semilogy(Eb_N0_dB,theoryBerAWGN,'om-','LineWidth',2); hold on
semilogy(Eb_N0_dB,theoryBerAWGN,'*g--','LineWidth',2); hold on
semilogy(Eb_N0_dB,nErr1,'ok--','LineWidth',2); hold on
semilogy(Eb_N0_dB,nErr2,'*r--','LineWidth',2); hold on
% semilogy(Eb_N0_dB,BER_Tx1Rx1_DIV,'ok-','LineWidth',2); hold on
% semilogy(Eb_N0_dB,BER_Tx2Rx2_DIV,'*m--','LineWidth',2); hold on
semilogy(Eb_N0_dBzf,theoryBer_nRx1,'pb-','LineWidth',2); hold on
semilogy(Eb_N0_dBzf,simBer,'pr--','LineWidth',2); hold on
semilogy(Eb_N0_dB,theoryBer,'cd--','LineWidth',2);
axis([0 20 10^-5 0.5])
grid on
legend('Rx1 proposed-ideal', 'Rx2 proposed-ideal','Rx1 proposed practical','Rx2 proposed practical','theory (nTx=1,nRx=1)','sim (nTx=2, nRx=2, ZF)','Rayleigh-Simulation','fontweight','bold');
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
% end













