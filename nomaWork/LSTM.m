%https://www.mathworks.com/help/deeplearning/ug/classify-sequence-data-using-lstm-networks.html

%Load data
[XTrain,YTrain] = japaneseVowelsTrainData;
XTrain(1:5)

%%
%Visualize the first time series in a plot. Each line corresponds to a feature.
figure
plot(XTrain{1}')
xlabel("Time Step")
title("Training Observation 1")
numFeatures = size(XTrain{1},1);
legend("Feature " + string(1:numFeatures),'Location','northeastoutside')

%% Prepare Data for Padding
% Get the sequence lengths for each observation.
numObservations = numel(XTrain);
for i=1:numObservations
    sequence = XTrain{i};
    sequenceLengths(i) = size(sequence,2);
end
% Sort the data by sequence length.
[sequenceLengths,idx] = sort(sequenceLengths);
XTrain = XTrain(idx);
YTrain = YTrain(idx);
% View the sorted sequence lengths in a bar chart.
%%
figure
bar(sequenceLengths)
ylim([0 30])
xlabel("Sequence")
ylabel("Length")
title("Sorted Data")

%% Define LSTM Network Architecture

inputSize = 12;
numHiddenUnits = 100;
numClasses = 9;

layers = [ ...
    sequenceInputLayer(inputSize)
    bilstmLayer(numHiddenUnits,'OutputMode','last')
    fullyConnectedLayer(numClasses)
    softmaxLayer
    classificationLayer]
%% Specify the training options
maxEpochs = 100;
miniBatchSize = 27;

options = trainingOptions('adam', ...
    'ExecutionEnvironment','cpu', ...
    'GradientThreshold',1, ...
    'MaxEpochs',maxEpochs, ...
    'MiniBatchSize',miniBatchSize, ...
    'SequenceLength','longest', ...
    'Shuffle','never', ...
    'Verbose',0, ...
    'Plots','training-progress');
%% Train LSTM Network
net = trainNetwork(XTrain,YTrain,layers,options);

%% Test LSTM Network
% Load test data
[XTest,YTest] = japaneseVowelsTestData;
XTest(1:3)

% Sort the test data by sequence length
numObservationsTest = numel(XTest);
for i=1:numObservationsTest
    sequence = XTest{i};
    sequenceLengthsTest(i) = size(sequence,2);
end
[sequenceLengthsTest,idx] = sort(sequenceLengthsTest);
XTest = XTest(idx);
YTest = YTest(idx);

% Classify the test data. 
miniBatchSize = 27;
YPred = classify(net,XTest, ...
    'MiniBatchSize',miniBatchSize, ...
    'SequenceLength','longest');
%%
acc = sum(YPred == YTest)./numel(YTest)














