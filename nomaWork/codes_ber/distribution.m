% 
% %https://www.mathworks.com/help/stats/histfit.html
% H1 = histfit(channel_H,100,'generalized extreme value'); hold on;
% set(H1(2),'color','m')
% % delete(H1(1))
% H1 = histfit(channel_H,100,'weibull'); hold on;
% set(H1(2),'color','g')
% delete(H1(1))
% H1 = histfit(channel_H,100,'nakagami'); hold on;
% set(H1(2),'color','b')
% delete(H1(1))
% H1 = histfit(channel_H,100,'gamma');
% delete(H1(1))
% set(H1(2),'color','r'); 
% ys = lines
% yticks([-1 -0.8 -0.2 0 0.2 0.8 1])


mu = mean(channel_H);
sigma = std(channel_H);


figure(1)
[f,x]=hist(channel_H,200);
y=f/trapz(x,f);
bar(x,y); hold on;

% x = 0:0.1:20;
% pd1 = makedist('Gamma','a',mu,'b',sigma);
% y = pdf(pd1, x);
% plot(x,y, 'k'); hold on;
% 
% pd2 = makedist('Weibull','a',mu,'b',sigma);
% y = pdf(pd2, x);
% plot(x,y, 'r'); hold on;
% 
% pd3 = makedist('Nakagami','mu',mu,'omega',sigma);
% y = pdf(pd3, x);
% plot(x,y, 'b'); hold on;
% 
% pd4 = makedist('GeneralizedExtremeValue','mu',mu,'sigma',sigma); 
% y = pdf(pd4, x);
% plot(x,y, 'g');
% ylabel('Probability Density');
% xlabel('Value');
% grid on;
% legend('Empirical','Gamma','Weibull','Nakagami','GeneralizedExtremeValue');
% 
% 

