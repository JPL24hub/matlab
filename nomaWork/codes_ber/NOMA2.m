clear all;
clc;
%%% NOMA parameters 
P = 1;
G1 = 10;
G2 = 10;
count = 1;
for alpha = 0:0.01:1 %power splitting factor
    P1 = P*alpha;
    P2 = P - P1;
    R1(count) = log2(1 + P1*G1);
    R2(count) = log2(1 + P2*G2/(P1*G2 + 1));
    count = count + 1;
end

hold on;
plot (R1,R2,'r');
grid on;
count = 1;

% 
% for alpha = 0:0.01:1 %bandwidth splitting factor 
%     P1 = P/2;
%     P2 = P/2;
%     R1(count) = alpha*log2(1 + P1*G1/alpha); R2(count) = (1-alpha)*log2(1 + P2*G2/(1-alpha)); count = count + 1;
% end
% hold on;
% plot(R1,R2,'k');
% xlabel('Rate of user 1 (bps/Hz)');
% ylabel('Rate of user 2 (bps/Hz)');
% grid on;
% box on;
% legend('NOMA','OFDMA')