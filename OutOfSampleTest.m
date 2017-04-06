%%
percentOut =  0.1;
percentIn = 1 - percentOut;
rng('shuffle');
clear ZOut ZIn;
for i = 1:nGroups
    % Sample size per class
    nSamples = size(Z{i},1);
    % Sample out size
    nSamplesOut = ceil(nSamples * percentOut);
    % Random select samples to stay out (test set)
    outSamples{i} = randperm(nSamples, nSamplesOut)';
    ZOut{i,1} = Z{i}(outSamples{i},:);
    % Selects samples that will be used to train and cross-validate the model
    samplesIndexes = (1:nSamples)';
    aux = ~ismember(samplesIndexes, outSamples{i});
    inSamples{i} = samplesIndexes(aux);
    ZIn{i,1} = Z{i}(inSamples{i},:);
end
%% Organize data
trainData = cell2mat(ZIn);
testData = cell2mat(ZOut);
% Generate labels, pick another subject data and put everything together
trainLabels = []; testLabels = [];
trainSubjectData = []; testSubjectData = [];
for i=1:nGroups
    trainSubjectData = [trainSubjectData ; FEATMATRIX(iia{i}(inSamples{i}), 1:5)];
    testSubjectData = [testSubjectData ; FEATMATRIX(iia{i}(outSamples{i}), 1:5)];
    
    trainLabels = [trainLabels ; FEATMATRIX(iia{i}(inSamples{i}),4) + ...
        FEATMATRIX(iia{i}(inSamples{i}),5) + 1];
    testLabels = [testLabels ; FEATMATRIX(iia{i}(outSamples{i}),4) + ...
        FEATMATRIX(iia{i}(outSamples{i}),5) + 1];
end
% clear x1;
% x1In = [s Z1InLabels];

%% Sammon
opts = sammon;
opts.TolFun = 1e-50;
opts.MaxHalves = 1000;
opts.MaxIter = 1000;
opts.LearningRate = 0.5;
[s, mappingInfo] = compute_mapping(trainData, 'Sammon', 2, opts);

figure;
scatter(wholeDataSammon(:,1), wholeDataSammon(:,2), 80, wholeLabelData);

% Train neural network
% Run PCA on high dimensional data and keep 90% of variance

% Simulate neural network
netOutSammon = sim(net_n35_mse28, wholeDataPCASammon')';
% Plot nn result
figure;
scatter(netOutSammon(:,1), netOutSammon(:,2), 80, wholeLabelData);
%% Parametric t-SNE

initial_dims = 30;
%perplexity = [30 30 30 30 30 30];
perplexity = [5 30 50];
max_iter = [1000 1000 1000];
%max_iter = [500];
learning_rate = 300;
% Network structure
layers = [15 15 60 2];
%layers = [500 500 2000 2];
diary on;
datetime('now')
perplexity
max_iter
for i=1:length(perplexity)
    tic
    % Train the parametric t-SNE network
    [network, err] = train_par_tsne(trainData, trainLabels, testData, testLabels, layers, 'CD1', ...
        max_iter(i), perplexity(i));
    % Construct training and test embeddings
    mappedTrain = run_data_through_network(network, trainData);
    mappedTest  = run_data_through_network(network, testData);
    % Compute 1-NN error and trustworthiness
    disp(['i: ' num2str(i) ' 1-NN error: ' num2str(knn_error(mappedTrain, trainLabels, mappedTest, testLabels, 1))]);
    disp(['i: ' num2str(i) ' Trustworthiness: ' num2str(trustworthiness(testData, mappedTest, 12))]);
    % Plot test embedding
    figure;
    scatter(mappedTest(:,1), mappedTest(:,2), 60, testLabels);
    title(['Embedding of test data - Perp: ' num2str(perplexity(i)) ' - MaxIter: ' num2str(max_iter(i))]);
    % Plot train embedding
    figure;
    scatter(mappedTrain(:,1), mappedTrain(:,2), 60, trainLabels);
    title(['Embedding of train data - Perp: ' num2str(perplexity(i)) ' - MaxIter: ' num2str(max_iter(i))]);
    toc
    disp(' ');
    disp('--------------------------------------------------------------');
    disp(' ');
end
diary off;
%% Default t-SNE
initial_dims = 30;
perplexity = 20;
max_iter = 1000;
learning_rate = 300;
[s, mappingInfo] = compute_mapping(trainData, 'tSNE', 2, initial_dims, perplexity, max_iter, learning_rate, trainLabels);

figure;
scatter(wholeDataTsne(:,1), wholeDataTsne(:,2), 80, wholeLabelData);

% Train neural network
% Run PCA on high dimensional data and keep 90% of variance
[wholeDataPCA, mappingInfoPCA] = compute_mapping(wholeData, 'PCA', 0.9);
mappingInfoPCA.error = -1;

% Simulate neural network
netOutTsne = sim(net_n30_mse4_tsne, wholeDataPCA')';
% Plot nn result
figure;
scatter(netOutTsne(:,1), netOutTsne(:,2), 80, wholeLabelData);

%% PCA

%% In samples plot...
Y = [];
for i=1:nGroups
    Y = [Y ; FEATMATRIX(iia{i}(inSamples{i}), 1:5)];
end
Y = [Y s];
% Scatter plot...
%%
testPoints = out_of_sample_est(testData, trainData, s);
%[testPoints, mappingInfo2] = compute_mapping(Z1Out, 'Sammon', 2, opts);
%% Out of samples plot...
Y = [];
for i=1:nGroups
    Y = [Y ; FEATMATRIX(iia{i}(outSamples{i}), 1:5)];
end
Y = [Y testPoints];
% Scatter plot...
