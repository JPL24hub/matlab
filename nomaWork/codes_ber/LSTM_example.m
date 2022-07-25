

N=1000;
A = randn(1, N);
B = randn(1, N);
C = randn(1, N);

for i=1:numel(A)
    Y(i) = (A(i)+B(i))/C(i);
end



%% Load the data
figure
plot(Y)
xlabel("Month")
ylabel("Cases")
title("Monthly Cases of Chickenpox")

%Split the data
numTimeStepsTrain = floor(0.9*numel(Y));

YTrain = Y(1:numTimeStepsTrain+1);
YTest = Y(numTimeStepsTrain+1:end);

ATrain = A(1:numTimeStepsTrain+1);
ATest = A(numTimeStepsTrain+1:end);

BTrain = B(1:numTimeStepsTrain+1);
BTest = B(numTimeStepsTrain+1:end);

CTrain = C(1:numTimeStepsTrain+1);
CTest = C(numTimeStepsTrain+1:end);

%% Standardize Data

myData = [YTrain ATrain BTrain CTrain];

mu = mean(myData);
sig = std(myData);


%% Prepare Predictors and Responses
ATrain = (ATrain - mu) / sig;
BTrain = (BTrain - mu) / sig;
CTrain = (CTrain - mu) / sig; %training sequences without the final time step

XTrain = [ATrain;BTrain;CTrain];
YTrain = (YTrain - mu) / sig;
%% Define LSTM Network Architecture

numFeatures = 3;
numResponses = 1;
numHiddenUnits = 200;

layers = [ ...
    sequenceInputLayer(numFeatures)
    lstmLayer(numHiddenUnits)
    fullyConnectedLayer(numResponses)
    regressionLayer];

% Specify the training options. Set the solver to 'adam' and train for 250 epochs. 
%To prevent the gradients from exploding, set the gradient threshold to 1. Specify 
% the initial learn rate 0.005, and drop the learn rate after 125 epochs by multiplying by a factor of 0.2.
options = trainingOptions('adam', ...
    'MaxEpochs',250, ...
    'GradientThreshold',1, ...
    'InitialLearnRate',0.005, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',125, ...
    'LearnRateDropFactor',0.2, ...
    'Verbose',0, ...
    'Plots','training-progress');
%% Train LSTM Network
net = trainNetwork(XTrain,YTrain,layers,options);

% %% Test
% % standadize first
% dataTestStandardized = (YTest - mu) / sig;
% XTest = dataTestStandardized(1:end-1);
% 
% net = predictAndUpdateState(net,XTrain);
% [net,YPred] = predictAndUpdateState(net,YTrain(end));
% 
% numTimeStepsTest = numel(XTest);
% for i = 2:numTimeStepsTest
%     [net,YPred(:,i)] = predictAndUpdateState(net,YPred(:,i-1),'ExecutionEnvironment','cpu');
% end
% 
% %Unstandardize the predictions using the parameters calculated earlier.
% YPred = sig*YPred + mu;
% 
% %Calculate the RMSE from the unstandardized predictions.
% YTest = YTest(2:end);
% rmse = sqrt(mean((YPred-YTest).^2))
% 
% %% Plots
% %Plot the training time series with the forecasted values.
% 
% figure
% plot(YTrain(1:end-1))
% hold on
% idx = numTimeStepsTrain:(numTimeStepsTrain+numTimeStepsTest);
% plot(idx,[Y(numTimeStepsTrain) YPred],'.-')
% hold off
% xlabel("Month")
% ylabel("Cases")
% title("Forecast")
% legend(["Observed" "Forecast"])
% 
% % Compare the forecasted values with the test data.
% figure
% subplot(2,1,1)
% plot(YTest)
% hold on
% plot(YPred,'.-')
% hold off
% legend(["Observed" "Forecast"])
% ylabel("Cases")
% title("Forecast")
% 
% subplot(2,1,2)
% stem(YPred - YTest)
% xlabel("Month")
% ylabel("Error")
% title("RMSE = " + rmse)
% 
% 
