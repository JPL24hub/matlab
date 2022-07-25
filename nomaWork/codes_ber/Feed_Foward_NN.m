

% inputData = cell2mat(channelData(1));
% outputData = cell2mat(RData(1));

N=64;
inputDataA = diag(randn(1, N));
inputDataB = diag(randn(1, N));
inputDataC = diag(randn(1, N));
inputDataD = diag(randn(1, N));
inputDataE = diag(randn(1, N));
inputDataF = diag(randn(1, N));



outputData1 = ((inputDataE*inputDataC)-(inputDataA*inputDataF)); %/((inputDataB*inputDataC)-(inputDataA*inputDataD));
outputData2 = ((inputDataB*inputDataC)-(inputDataA*inputDataD));
outputData3 = (inputDataF-(inputDataA*inputDataC));


%% net3
outputData31 = diag(outputData3)';
inputData31 = [diag(inputDataA)';diag(inputDataC)';diag(inputDataF)'];

% [nomInput, nomOutput, sig, mu] = normalize_unormalize(inputData31, outputData31, 1, 0, 0);

%%

net3 = feedforwardnet([5 10 5]);
net3.trainFcn = 'trainbr';
[net3, tr3] = train(net3, inputData31, outputData31);
view(net3)

%%
[Input, Output, sig, mu] = normalize_unormalize(nomInput, nomOutput, 2, sig, mu);
 figure(1)
 plot(outputData31)
 figure(2)
 plot(nomOutput)
 figure(3)
 plot(Output)






















%% net 1
outputData1 = diag(outputData1)';
inputData1 = [diag(inputDataA)';diag(inputDataC)';diag(inputDataE)';diag(inputDataF)'];

net1 = feedforwardnet([5 10 5]);
net1.trainFcn = 'trainbr';
[net1, tr1] = train(net1, inputData1, outputData1);
view(net1)

%% net 2
outputData2 = diag(outputData2)';
inputData2 = [diag(inputDataA)';diag(inputDataB)';diag(inputDataC)';diag(inputDataD)'];

net2 = feedforwardnet([5 10 5]);
net2.trainFcn = 'trainbr';
[net2, tr2] = train(net2, inputData2, outputData2);
view(net2)

 



%% important links
%https://www.mathworks.com/matlabcentral/answers/219775-multiple-input-feedforwardnet-neural-network-toolbox