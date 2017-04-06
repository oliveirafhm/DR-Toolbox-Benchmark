function [net, tr, e, performance, r] = OutOfSampleNN(x, t, trainFcn, ...
    hiddenLayerSize, trainRatio, valRatio, testRatio)
%OUTOFSAMPLENN Creates and trains a neural network to learn embedding
%(low-dimensional map)
%   Date......: March 09, 2017 
%   Author....: F?bio Henrique M Oliveira (oliveirafhm@gmail.com)
%   Parameters:
%               x -> Training sample data
%               t -> NN Targets
%               trainFcn -> Network training function (see matlab docs)
%               hiddenLayerSize -> Quantity of neurons in hidden layer
%               trainRatio/valRatio/testRatio -> percent of data to be used
%               in each NN training step
%   Return....:
%               net -> Trained network
%               tr -> Training record (epoch and perf)
%               e -> Errors
%               perfomance -> Overall and test set performance (mse)
%               r -> Correlation coefficient (R-value) between the outputs
%               and targets (overall and test set)

% Create a Fitting Network
net = fitnet(hiddenLayerSize,trainFcn);
% Test
net.trainParam.epochs = 1000;
% Setup Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = trainRatio/100;
net.divideParam.valRatio = valRatio/100;
net.divideParam.testRatio = testRatio/100;
% Training parameters
%net.trainParam.showWindow = false;
% Train the Network (TODO: test use gpu)
[net,tr] = train(net,x',t');%,'useParallel','yes');
% Test the Network
y = net(x');
e = gsubtract(t',y);
performance.overall = perform(net,t',y);
r.overall = regression(t',y,'one');
% Performance on the test set
tInd = tr.testInd;
y = net(x(tInd,:)');
performance.test = perform(net,t(tInd,:)',y);
r.test = regression(t(tInd,:)',y,'one');
end

