

clear all

Eb_N0_dB = 0:1:40; % multiple Eb/N0 values
n_iter = 1000;
symbol_len=64;



for ii = 1:length(Eb_N0_dB)
    BER1=0;
    bittot=0;
    
    for k = 1:n_iter
        
        x1 = rand(symbol_len,1)>0.5; % generating 0,1 with equal probability
        s1 = 2*x1-1; % BPSK modulation 0 -> -1; 1 -> 0
        
        n1 = 1/sqrt(2)*[randn(1,symbol_len) + 1i*randn(1,symbol_len)];     
        h11 = 1/sqrt(2)*[randn(1,symbol_len) + 1i*randn(1,symbol_len)]; % Rayleigh channel
        
        % receiver - hard decision decoding
        y1 = diag(h11)*diag(s1) + diag(n1)*10^(-Eb_N0_dB(ii)/20);
        y1hat = y1/diag(h11);
        
        x1Hat = real(diag(y1hat))>0;
    
        
        BER1 = BER1 + size(find([x1 - x1Hat]),1);
        bittot = bittot + symbol_len;
        
    end
    nErr1(ii) = BER1/bittot;
end
%%
semilogy(Eb_N0_dB,nErr1,'bp-','LineWidth',2);hold on
grid on
legend('TX1');
xlabel('Average Eb/No,dB');
ylabel('Bit Error Rate');
% end













