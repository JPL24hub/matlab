
clc, clear all, warning off;

x = randi([0,1],1000,1);
%Mapping
for i=1:length(x)
    if x(i)==1
        xmode(i)=1;
    else
        xmode(i)=-1;
    end
end
% scatterplot(xmode) %observe constalation points
BER_sim =[]; %simulated BER
% scatterplot(y)
dt=2
for snr=0:dt:20  %in db
    y = awgn(complex(xmode),snr); %modulated sysmbols transmited through awgm
    for i=1:length(y)
        if y(i)>=0
            detValue(i,1)=1;%detected value
        else
            detValue(i,1)=0;
        end
    end
[nor,ber] = biterr(x,detValue);
BER_sim= [BER_sim ber];
end

snr=0:dt:20;

figure(1);
semilogy(snr,BER_sim,'^r-','LineWidth',2); 
xlabel('SNR (dB)');
ylabel('Bit Error Rate(BER)');
legend('BER');




