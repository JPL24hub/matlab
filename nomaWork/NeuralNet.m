a1=[1 2 3 4 5 6 7 8 9 10]'; 
b1=[3 4 5 6 7 8 9 10 11 1]'; 
c1=[1 4 2 5 3 6 4 8 9 10]'; 
d1=[2 4 5 6 7 9 8 10 1 11]'; 
e1=[1 6 3 5 4 7 8 9 10 11]'; 
f1=[10 9 8 7 1 2 3 4 5 6 ]';
o1=[50 51 52 54 58 60 65 66 69 70]';

o2=[50 51 52 54 58 60 65 66 69 70]';

o = [o1 o2];

% Create input matrix and output vector 
inputs=[a1 b1 c1 d1 e1 f1]'; targets=o';
%%
inputs = [a b c d x_1' x_2']';
o = [R1 R2];
targets=o';
%%
% Create a Fitting Network 
hiddenLayerSize = 4; 
net = fitnet(hiddenLayerSize);

% Choose Input and Output Pre/Post-Processing Functions % For a list of all processing functions type: help nnprocess 
net.inputs{1}.processFcns = {'mapminmax'}; 
net.outputs{2}.processFcns = {'mapminmax'};

% Setup Division of Data for Training, Validation, Testing % For a list of all data division functions type: help nndivide 
net.divideFcn = 'dividerand'; % Divide data randomly 
net.divideMode = 'sample'; % Divide up every sample 
net.divideParam.trainRatio = 70/100; 
net.divideParam.valRatio = 15/100; 
net.divideParam.testRatio = 15/100;

net.trainParam.epochs=100;

% For help on training function 'trainlm' type: help trainlm % For a list of all training functions type: help nntrain 
net.trainFcn = 'trainlm'; % Levenberg-Marquardt

% Choose a Performance Function % For a list of all performance functions type: help nnperformance 
net.performFcn = 'mse'; % Mean squared error
% Choose Plot Functions % For a list of all plot functions type: help nnplot 
net.plotFcns = {'plotperform','plottrainstate','ploterrhist', 'plotregression', 'plotfit'};
% Train the Network 
[net,tr] = train(net,inputs,targets);

% Test the Network 
outputs = net(inputs); 
errors = gsubtract(targets,outputs); 
performance = perform(net,targets,outputs)

% Recalculate Training, Validation and Test Performance 
trainTargets = targets .* tr.trainMask{1}; 
valTargets = targets .* tr.valMask{1}; 
testTargets = targets .* tr.testMask{1}; 
trainPerformance = perform(net,trainTargets,outputs) 
valPerformance = perform(net,valTargets,outputs) 
testPerformance = perform(net,testTargets,outputs)

% View the Network view(net)
%figure, plotfit(net,inputs,targets)