
% function[net1,net2,net3] = WP2_FFW(channelData, RData)
    inputData = [];
    outputData = [];
   for i=1:700 
    inputData = [inputData cell2mat(channelData(i))];
    outputData = [outputData cell2mat(RData(i))];
   end

%% net 1
    
    inputData1 = [inputData(5,:);inputData(3,:);inputData(1,:);inputData(6,:)];
    outputData1 = outputData(3,:);
    
    net1 = feedforwardnet([5 10 5]);
    net1.trainFcn = 'trainbr';
    [net1, tr1] = train(net1, inputData1, outputData1);
    view(net1)
%     %% test
%     inputData = cell2mat(channelData(i+1));
%     inputData1 = [inputData(5,:);inputData(3,:);inputData(1,:);inputData(6,:)];
%     outputData = cell2mat(RData(i+1));
%     Pred1 = net1(inputData1);
%   
%     plot(Pred1,'r')
%     hold on;
%     plot(outputData(3,:))
%% net 2

    inputData2 = [inputData(2,:);inputData(3,:);inputData(1,:);inputData(4,:)];
    outputData2 = outputData(4,:);

    net2 = feedforwardnet([5 10 5]);
    net2.trainFcn = 'trainbr';
    [net2, tr2] = train(net2, inputData2, outputData2);
    view(net2)

%% net3 R1 NUME
    inputData3 = inputData;
    outputData3 = outputData(5,:);

    net3 = feedforwardnet([5 10 5]);
    net3.trainFcn = 'trainbr';
    [net3, tr3] = train(net3, inputData3, outputData3);
    view(net3)

%% net4 R1 DENOM
    inputData4 = [inputData(2,:);inputData(3,:);inputData(1,:);inputData(4,:)];
    outputData4 = outputData(6,:);

    net4 = feedforwardnet([5 10 5]);
    net4.trainFcn = 'trainbr';
    [net4, tr4] = train(net4, inputData4, outputData4);
    view(net3)
    



%% important links
%https://www.mathworks.com/matlabcentral/answers/219775-multiple-input-feedforwardnet-neural-network-toolbox