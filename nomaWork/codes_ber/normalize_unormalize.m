
% function [normA] = normalize_unormalize(input)

% if signal == 1
%     alldata = reshape(input,[1,numel(input)]);
%     size(alldata)
%     mu = mean(alldata)
%     sig = std(alldata)
%     
%     nomInput =(input - mu) / sig;
    
    normA = input - min(input(:));
    normA = normA ./ max(normA(:));

%     nomOutput = (output - mu) / sig;
% % else
%     nomInput = sig*input + mu;
%     nomOutput = (output - mu) / sig;
% end
%%
berr1 = BERu11;
ber1r2 = BERu21;
berrr2 = BERIu11;

% end


