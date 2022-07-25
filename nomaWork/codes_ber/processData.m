
function[normalizeInput, normalizeIOutput, sig, mu] = processData(channelData,RData)

% channelData = channelData';
len1 = length(channelData);

% RData =  RData';

channelDatamat = cell2mat(channelData); %channel data to matrix
RDatamat = cell2mat(RData); %channel data to matrix


Data = [channelDatamat;RDatamat]; %combine all data for normalization


channelDataArr = reshape(Data,[],1); %make new array from datato get mu and sig
mu = mean(channelDataArr);
sig = std(channelDataArr);



nomalizedChan = (channelDatamat - mu) / sig; %normalize data
nomalizedR = (RDatamat - mu) / sig; %normalize data



rowDist1 = ones(1,(size(nomalizedChan,1)/size(channelData{1}, 1)))*size(channelData{1}, 1);
normalizeInput = mat2cell(nomalizedChan,rowDist1); %6000 by 1 cell 

rowDist2 = ones(1,(size(nomalizedR,1)/size(RData{1}, 1)))*size(RData{1}, 1);
normalizeIOutput = mat2cell(nomalizedR,rowDist2); %6000 by 1 cell 



% const = 5;
% j=1;
% normalizeInput = {}; 
% for i=1:len1 %take the first 3600 by 64 to get norm input to cell
%     normalizeInput{i} = nomalizedChan(j:j+const,:);
%     j=j+6;
% end
% 
% const = 1;
% j=1;
% normalizeIOutput = {};
% for i=1:len1 %take the last from j to last to get norm input to cell
%     normalizeIOutput{i} = nomalizedR(j:j+const,:);  %
%     j=j+2;  
% end
% 
% normalizeInput = normalizeInput';
% normalizeIOutput = normalizeIOutput';
end
