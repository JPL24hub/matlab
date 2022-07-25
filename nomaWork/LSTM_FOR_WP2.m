
 function [newNet, i] = LSTM_FOR_WP2(norm_channelData, norm_RData)
%Load data
XTrain = norm_channelData(1);
YTrain = norm_RData(1);

%https://www.mathworks.com/help/deeplearning/ref/nnet.cnn.layer.lstmlayer.html

%% Define LSTM Network Architecture
numFeatures = 6;
numResponses = 2;
numHiddenUnits1 = 384;
numHiddenUnits2 = 192;
numHiddenUnits3 = 128;
numHiddenUnits4 = 64;

layers = [ ...
    sequenceInputLayer(numFeatures)
    lstmLayer(numHiddenUnits1)
    dropoutLayer(0.8)
    
    lstmLayer(numHiddenUnits2)
    dropoutLayer(0.8)
    
    lstmLayer(numHiddenUnits3)
    dropoutLayer(0.8)
    
    lstmLayer(numHiddenUnits4)
    dropoutLayer(0.8)
    
    fullyConnectedLayer(numResponses)
    regressionLayer];

%% Specify the training options
maxEpochs = 2;
miniBatchSize = 48;

options = trainingOptions('adam', ...
    'ExecutionEnvironment','auto', ...
    'MaxEpochs',maxEpochs, ...
    'MiniBatchSize',miniBatchSize, ...
    'Verbose',false);


%% Train LSTM Network

net = trainNetwork(XTrain,YTrain,layers,options); 

% Extract layers from the trained network
newLayers = net.Layers;

len = 1000;%length(norm_channelData)-1

    for i = 2 : len
  
        fprintf('%d out of %d\n', i, len)
        XTrain = norm_channelData(i);
        YTrain = norm_RData(1);
    % retrain network using layers from previous network
        newNet = trainNetwork(XTrain, YTrain, newLayers, options);
        newLayers = newNet.Layers;
    end
    
    i
    
    options = trainingOptions('adam', ...
    'ExecutionEnvironment','auto', ...
    'MaxEpochs',maxEpochs, ...
    'MiniBatchSize',miniBatchSize, ...
    'Verbose',false, ...
    'Plots','training-progress');

    XTrain = norm_channelData(i);
    YTrain = norm_RData(i);
    
    newNet = trainNetwork(XTrain, YTrain, newLayers, options);
    
    i=i+1;
    
    


 end



