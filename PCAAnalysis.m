% By F?bio Henrique (oliveirafhm@gmail.com)
% Run before Classification.m
% Dimensional reduction (t-SNE, Sammon and PCA)
% Taking in account each task for all groups
% 04/2017
%% Global vars
if useConfigData == 0
    saveData = inputdlg('Dou you want to save all data related with DR and Classification? y/n');
    saveData = saveData{1};
elseif useConfigData == 1
    saveData = expConfigData(mainCounter).SaveData;
end

%%
% Sinal filtrado, normalizado e devidamente centralizado em torno do zero
load 'FeatMatrix9.mat';
% Frequencia instantanea  (Baseada na transforma de Hilbert)
load 'FeatMatrix10.mat';
% Amplitude instantanea   (Baseada na transforma de Hilbert)
load 'FeatMatrix11.mat';

%           1    2    3    4       5       6       7
dataSet = {'FS','IF','IA','FS-IF','FS-IA','IF-IA','FS-IF-IA'};
if useConfigData == 0
    dsChoice = inputdlg('Which dataset do you want to use? 1-7 (check dataSet var)');
    dsChoice = str2double(dsChoice{1});
elseif useConfigData == 1
    dsChoice = expConfigData(mainCounter).DataSet;
end
selectedDataSet = dataSet{dsChoice};

switch selectedDataSet
    case dataSet{1}
        FEATMATRIX = FeatMatrix9;
    case dataSet{2}
        FEATMATRIX = FeatMatrix10;
    case dataSet{3}
        FEATMATRIX = FeatMatrix11;
    case dataSet{4}
        FEATMATRIX = [FeatMatrix9 FeatMatrix10(:,6:end)];
    case dataSet{5}
        FEATMATRIX = [FeatMatrix9 FeatMatrix11(:,6:end)];
    case dataSet{6}
        FEATMATRIX = [FeatMatrix10 FeatMatrix11(:,6:end)];
    case dataSet{7}
        FEATMATRIX = [FeatMatrix9 FeatMatrix10(:,6:end) FeatMatrix11(:,6:end)];
    otherwise
        error('Unknown data set arrangement.');
end

%%
% 1 = non parkinson | 2 = parkinson | 3 = all dbs |
% 4 = x dbs intensity (not implemented)

% Norm each group separately (0)
% Normalize all groups using mean and std from non parkinson group (1)
% Norm all groups together (2)
normMethod = 2;
nGroups = 3;
% Choose which task will be analyzed (0 to all)
% 1 - Finger taps | 2 - Finger to nose | 3 - Supination and pronation
% 4 - Rest
if useConfigData == 0
    task = inputdlg('Type task number? 1 / 2 / 3 / 4');
    task = str2double(task{1});
elseif useConfigData == 1
    task =  expConfigData(mainCounter).TaskNumber;
end
% Before normalization
z = cell(nGroups,1);
% After normalization
Z = cell(nGroups,1);
clear t;
for i=1:nGroups
    if i == 1
        % Pick only non parkinson subjects
        if task == 0
            iia{i} = find(FEATMATRIX(:,4) == 0);
        else
            iia{i} = find(FEATMATRIX(:,4) == 0 & FEATMATRIX(:,2) == task);
        end
        t{i} = 'S_{H}';
        z{i} = FEATMATRIX(iia{i},6:end);
        Z{i} = zeros(size(z{i}));
    elseif i == 2
        % Pick only parkinson subjects
        if task == 0
            iia{i} = find(FEATMATRIX(:,4) == 1 & FEATMATRIX(:,5) == 0);
        else
            iia{i} = find(FEATMATRIX(:,4) == 1 & FEATMATRIX(:,5) == 0 ...
                & FEATMATRIX(:,2) == task);
        end
        t{i} = 'S_{PD}';
        z{i} = FEATMATRIX(iia{i},6:end);
        Z{i} = zeros(size(z{i}));
    elseif i == 3
        % Pick only dbs subjects
        if task == 0
            iia{i} = find(FEATMATRIX(:,4) == 1 & FEATMATRIX(:,5) == 1 );
        else
            iia{i} = find(FEATMATRIX(:,4) == 1 & FEATMATRIX(:,5) == 1 ...
                & FEATMATRIX(:,2) == task);
        end
        t{i} = 'S_{DBS}';
        z{i} = FEATMATRIX(iia{i},6:end);
        Z{i} = zeros(size(z{i}));
    end
    if normMethod == 0
        Z{i} = zscore(z{i});
    elseif normMethod == 1
        % Data normalization using mean and std from non parkinson group
        if i == 1
            [Z{i}, mu, sigma] = zscore(z{i});
        else
            for j=1:length(mu)
                Z{i}(:,j) = (z{i}(:,j) - mu(j)) / sigma(j);
            end
        end
    elseif normMethod == 2 && i == nGroups
        Z0 = cell2mat(z);
        Z0 = zscore(Z0);
        zAux = [0];
        for j=1:nGroups
            zAux = [zAux ; zAux(j)+length(iia{j})];
            Z{j} =  Z0(zAux(j)+1:zAux(j+1),:);
        end
        clear Z0;
    end
end
%% Split data (train and test sets)
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
% Generate labels, pick subject data and put everything together
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

%% PCA, Sammon and t-SNE projection
%addpath(genpath('Third party codes'));

% Group groups before run 2D projection?
group = 1;

computeMapAlgLabels = {'PCA','Sammon','t-SNE'};
if useConfigData == 0
    computeMapAlg = inputdlg('Choose DR algorithm? 1(pca) / 2(sammon) / 3(t-sne)');
    computeMapAlg = str2double(computeMapAlg{1});
elseif useConfigData == 1
    computeMapAlg = expConfigData(mainCounter).DRAlg;
end

% 2D data projection
% Group
if group == 1
    %     Y = [];
    clear x1 x2;
    Z1 = trainData;
    if computeMapAlg == 2
        opts = sammon;
        opts.TolFun = 1e-50;
        opts.MaxHalves = 1000;
        if useConfigData == 0
            opts.MaxIter = inputdlg('Type MaxIter param: (int between 1000 - 5000)');
            opts.MaxIter = str2double(opts.MaxIter{1});
            opts.LearningRate = inputdlg('Type LearningRate param: (between 0.3 - 0.7)');
            opts.LearningRate = str2double(opts.LearningRate{1});
        elseif useConfigData == 1
            opts.MaxIter = expConfigData(mainCounter).MaxIter;
            opts.LearningRate = expConfigData(mainCounter).LearningRate;
            if opts.MaxIter <= 0 || opts.LearningRate <= 0
                error('Invalid parameters (MaxIter or LearningRate) for Sammon algorithm.');
            end
        end
        [s, mappingInfo] = compute_mapping(Z1, 'Sammon', 2, opts);
    elseif computeMapAlg == 3
        % t-SNE projection algorithm parameters
        initial_dims = 30;
        % https://lvdmaaten.github.io/tsne/
        if useConfigData == 0
            perplexity = inputdlg('Type perplexity param: (int between 5 - 50)');
            perplexity = str2double(perplexity{1});
            max_iter = inputdlg('Type max_inter param: (int between 1000 - 5000)');
            max_iter = str2double(max_iter{1});
            learning_rate = inputdlg('Type learning_rate param: (between 300 - 700)');
            learning_rate = str2double(learning_rate{1});
        elseif useConfigData == 1
            perplexity = expConfigData(mainCounter).Perplexity;
            max_iter = expConfigData(mainCounter).MaxIter;
            learning_rate = expConfigData(mainCounter).LearningRate;
            if max_iter <= 0 || learning_rate <= 0 || perplexity <= 0
                error('Invalid parameters (MaxIter, LearningRate or Perplexity) for t-SNE algorithm.');
            end
        end
        [s, mappingInfo] = compute_mapping(Z1, 'tSNE', 2, initial_dims, perplexity, max_iter, learning_rate,[]);
        mappingInfo.iter = -1;
    elseif computeMapAlg == 1
        [s, mappingInfo] = compute_mapping(Z1, 'PCA', 2);
        mappingInfo.error = -1;
        mappingInfo.iter = -1;
    end
    % Reassembly of data array with each individual data + projections
    % instead of features
    % Train subject data + DR projection
    Y = [trainSubjectData s];
    % Projection + class (group)
    x1 = [s trainLabels];
    
    % Do not group (do the individual projection for each group)
    % Removed ...discontinued.
end
% Save NIter (Sammon only) and DRError in ExperimentsConfig xls file
if useConfigData == 1
    dRErrorCol = colHeadingsLetters{ismember(colHeadingsTxt,'DRError')};
    xlwrite(configFilePath, mappingInfo.error, 1, [dRErrorCol int2str(mainCounter+1)]);
    nIterCol = colHeadingsLetters{ismember(colHeadingsTxt,'NIter')};
    xlwrite(configFilePath, mappingInfo.iter, 1, [nIterCol int2str(mainCounter+1)]);
end

%% Define baseName, resultsXlsName and titleBaseName vars
taskNames = {'Finger taps','Finger to nose','Supination and pronation',...
    'Rest'};
% Var also used in Classification script
baseName = ['T' int2str(task) '_' taskNames{task} '_' ...
    computeMapAlgLabels{computeMapAlg} '_' selectedDataSet '_' date];
% Sammon
if computeMapAlg == 2
    resultsXlsName = ['_MI_' int2str(opts.MaxIter) ...
        '_LR_' num2str(opts.LearningRate)];
    baseName = [baseName resultsXlsName];
    % t-SNE
elseif computeMapAlg == 3
    resultsXlsName = ['_MI_' int2str(max_iter) ...
        '_LR_' num2str(learning_rate) '_PE_' ...
        int2str(perplexity)];
    baseName = [baseName resultsXlsName];
    %elseif computeMapAlg == 1
    % PCA does not have any extra parameters for the moment
end
titleBaseName = strrep(baseName, '_', ' ');
% Save BaseFileName in ExperimentsConfig xls file
if useConfigData == 1
    baseFileNameCol = colHeadingsLetters{ismember(colHeadingsTxt,'BaseFileName')};
    xlwrite(configFilePath, {baseName}, 1, [baseFileNameCol int2str(mainCounter+1)]);
end

%% Work folder
if useConfigData == 0
    workFolder = inputdlg('Type work folder name:');
    workFolder = workFolder{1};
elseif useConfigData == 1
    workFolder = expConfigData(mainCounter).WorkFolder;
end

addpath(genpath(workFolder));
startPath = pwd;
cd(workFolder);
d1 = dir([int2str(dsChoice) '*']);
cd(d1.name);
d2 = dir(['*' int2str(task)]);
cd(d2.name);
iterationPath = [workFolder '/' d1.name '/' d2.name];
cd(startPath);

%% Train and test NN
% Run PCA on high dimensional train data and keep 30 (or 90% of variance = x components) 
% first components (if DR method is not PCA)
if computeMapAlg ~= 1
%     [trainDataPCA, mappingInfoPCA] = compute_mapping(trainData, 'PCA', 30);
    [trainDataPCA, mappingInfoPCA] = compute_mapping(trainData, 'PCA', 0.9);
    % Map out of sample data (testData) to be used in the end (by SVM classifier)
    testDataPCA = out_of_sample(testData, mappingInfoPCA);
    % Bayesian Regularization backpropagation training function
    trainFcn = 'trainbr';
    % Hidden layer size
    if computeMapAlg == 2 % Sammon
        hiddenLayerSize = floor(mappingInfoPCA.no_dims * 0.5);
    elseif computeMapAlg == 3 % t-SNE
        hiddenLayerSize = floor(mappingInfoPCA.no_dims * 0.7);
    end
    % Setup Division of Data for Training (%), Validation (%), Testing (%)
    trainRatio = 75; % actually 90%
    valRatio = 15; % not used because of used training function
    testRatio = 10;
    % Train the Network
    nnTrainTrials = 5;
    bestNN = 1;
    for i=1:nnTrainTrials
        i
        [trainedN(i).net,trainedN(i).tr,trainedN(i).e,trainedN(i).performance,...
            trainedN(i).r] = OutOfSampleNN(trainDataPCA,s,trainFcn,...
            hiddenLayerSize,trainRatio,valRatio,testRatio);
        if trainedN(i).performance.test < trainedN(bestNN).performance.test
            bestNN = i;
        end
    end
% If DR method is PCA, we does not need train a model to out-of-sample
% process
elseif computeMapAlg == 1
    % Map out of sample data (testData) to be used in the end (by SVM classifier)
    mappingInfoPCA.no_dims = -1;
    testDataPCA = out_of_sample(testData, mappingInfo);
end
% Out-of-sample
if computeMapAlg ~= 1
    % Reduce dimensionality of new data using trained model (that learned
    % Sammon's or t-SNE function)
    oos.projection = trainedN(bestNN).net(testDataPCA')';
%     oos.performance = perform(trainedN(bestNN).net,??',oos.projection');
%     oos.r = regression(??',oos.projection','one');
elseif computeMapAlg == 1
    oos.projection = testDataPCA;
end
% Save trained NN info and stats, and related objects
if useConfigData == 1
    npcaCol = colHeadingsLetters{ismember(colHeadingsTxt,'NPCA')};
    xlwrite(configFilePath, mappingInfoPCA.no_dims, 1,...
            [npcaCol int2str(mainCounter+1)]);
    mseCol = colHeadingsLetters{ismember(colHeadingsTxt,'MSEOOS')};
    rCol = colHeadingsLetters{ismember(colHeadingsTxt,'ROOS')};
    if computeMapAlg ~= 1 % Sammon and t-SNE
        xlwrite(configFilePath, trainedN(bestNN).performance.test, 1,...
            [mseCol int2str(mainCounter+1)]);
        xlwrite(configFilePath, trainedN(bestNN).r.test, 1,...
            [rCol int2str(mainCounter+1)]);
    elseif computeMapAlg == 1 % PCA
        xlwrite(configFilePath, -1, 1, [mseCol int2str(mainCounter+1)]);
        xlwrite(configFilePath, -1, 1, [rCol int2str(mainCounter+1)]);
    end
end
if saveData == 'y'
    oosModelFileName = [baseName '_OosModelData.mat'];
    if computeMapAlg ~= 1
        bestTrainedN = trainedN(bestNN);
        save([iterationPath '/' oosModelFileName],'percentOut','trainData',...
            'testData','trainLabels','testLabels','trainSubjectData',...
            'testSubjectData','bestTrainedN','oos');
    elseif computeMapAlg == 1
        save([iterationPath '/' oosModelFileName],'percentOut','trainData',...
            'testData','trainLabels','testLabels','trainSubjectData',...
            'testSubjectData','oos');
    end
end

%% Plot DR method projection
baseFontSize = 17;
markerSize = 10;
figure;
hold on;
% Last color will be used by test set points
gColors = ['b','k','r','m'];
gMarkers = ['^','*','o'];
for i=1:nGroups
    if i == 1
        % Pick only non parkinson subjects
        Y1 = Y(Y(:,4) == 0,:);
    elseif i == 2
        % Pick only parkinson subjects
        Y1 = Y(Y(:,4) == 1 & Y(:,5) == 0,:);
    elseif i == 3
        % Pick only dbs subjects
        Y1 = Y(Y(:,4) == 1 & Y(:,5) == 1,:);
    end
    
    % Task x (var task was declared in the first section of this script)
    iib = find(Y1(:,2) == task);
    plot(Y1(iib,6), Y1(iib,7), [gColors(i) gMarkers(i)], ...
        'MarkerSize',markerSize);
end
% For to plot test set points
for i=1:nGroups
    iib = find(testLabels == i);
    plot(oos.projection(iib,1),oos.projection(iib,2), ...
        [gColors(end) gMarkers(i)], 'MarkerSize',markerSize);
end
title(titleBaseName, 'FontSize', baseFontSize+2);

legend(t{1},t{2},t{3},[t{1} '(TP)'],[t{2} '(TP)'],[t{3} '(TP)'],...
    'Location','Best','fontsize',baseFontSize);
xlabel('Dimension 1','fontsize',baseFontSize);
ylabel('Dimension 2','fontsize',baseFontSize);
set(gca,'fontsize',baseFontSize);
grid on;
hold off;

% Save figure
if saveData == 'y'
    figureName = [baseName '_ScatterPlot.fig'];
    savefig([iterationPath '/' figureName]);
    close;
end
