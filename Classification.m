% Test and train set classification stats
% Should be runned after PCAAnalysis
% By F?bio Henrique (oliveirafhm@gmail.com)
% Support Vector Machine
% 04/2017
%% Creates class variables
if useConfigData == 0
    genNewModel = inputdlg('Do you want to generate a new model? y/n');
    genNewModel = genNewModel{1};
elseif useConfigData == 1
    genNewModel = expConfigData(mainCounter).GenNewModel;
end

if genNewModel == 'y'
    x = x1;
elseif genNewModel == 'n'
    % Alert, you need to manual load model in this section
    [modelFileName, modelPathName] = uigetfile('.mat', ...
        'Select model file');
    modelFilePath = strcat(modelPathName, modelFileName);
    load(modelFilePath);
    xs = CVMdl.X;
    ys = CVMdl.Y;
    task = inputdlg('Type task number used to train this model? 1/2/3/4');
    task = str2double(task{1});
    ys = cellfun(@(x) str2double(x), ys,'UniformOutput', false);
    ys = cell2mat(ys);
    x = [xs ys];
end
if task == 1
    ts = 'Finger taps';
elseif task == 2
    ts = 'Finger to nose';
elseif task == 3
    ts = 'Supination and pronation';
elseif task == 4
    ts = 'Rest';
end

groupsNames = {'S_{H}','S_{PD}',...
    'S_{DBS}'};
% groupsNames = {'Non parkinson subjects','Parkinson subjects',...
%     'DBS subjects'};
classNumbers = unique(x(:,3));
nClass = length(classNumbers);
targetClassS = cell(length(x(:,3)),1); % targetClassSequence
for i = 1:length(targetClassS)
    targetClassS{i} = num2str(x(i,3));
end
classNames = unique(targetClassS);

%% Creates and trains SVM model
% TODO: Plot model boundaries (partially did, look DecisionBoundary.m)
if genNewModel == 'y'
    if useConfigData == 0
        % Default setup
        kernelFunction = 'gaussian';
        kernelScale = 0.35;
        codingDesign = 'onevsall';
    elseif useConfigData == 1
        kernelFunction = expConfigData(mainCounter).KernelFunction;
        kernelScale = expConfigData(mainCounter).KernelScale;
        cde = expConfigData(mainCounter).CodingDesign;
        if strcmpi(cde, 'ovo')
            codingDesign = 'onevsone';
        elseif strcmpi(cde, 'ova')
            codingDesign = 'onevsall';
        end
    end
    
    % Kernel chosen based on several tests done using classificationLearner
    % matlab tool.
    template = templateSVM(...
        'KernelFunction', kernelFunction, ...
        'PolynomialOrder', [], ...
        'KernelScale', kernelScale, ...
        'BoxConstraint', 1, ...
        'Standardize', true);
    
    % Init parallel pool if closed
    %     if isempty(gcp('nocreate')) == 1
    %         pool = parpool;
    %     end
    parOptions = statset('UseParallel',1);
    % TODO: Separate model creation and training from cross-validation
    CVMdl = fitcecoc(x(:,1:2),targetClassS, 'Learners', template, ...
        'FitPosterior',1, ...
        'ClassNames', classNames, 'Leaveout', 'on', 'Verbose', 1,...
        'Coding', codingDesign, 'Options', parOptions);
    
    % Only works in Matlab 2016b
    %     hpoo = struct('AcquisitionFunctionName',...
    %         'expected-improvement-plus');
    % 'HyperparameterOptimizationOptions', hpoo,
    %     CVMdl = fitcecoc(x(:,1:2),targetClassS, 'OptimizeHyperparameters','all',...
    %         'Options', parOptions);
    
    if saveData == 'y'
        modelFileName = [baseName '_LooSvmModel.mat'];
        save([iterationPath '/' modelFileName],'CVMdl','computeMapAlgLabels',...
            'computeMapAlg','dataSet','dsChoice','selectedDataSet','task');
    end
end
%% Statistical measures of the performance of a binary classification test
nMeasures = 9;
statisticalM = zeros(nClass, nMeasures);
overallSuccessRate = 0;
classificationStats = zeros(nClass,nClass);
rocMeasures = 5;
rocData = cell(nClass, rocMeasures);
% meanPScore = zeros(length(x), nClass);
% Predict responses for observations not used for training (CV model)
[label,~,PScore,Posterior] = kfoldPredict(CVMdl, 'Options', parOptions);

% Overall success rate for all classes (%)
overallTP = strcmp(targetClassS, label);
overallSuccessRate = 100*sum(overallTP)/length(overallTP)...
    + overallSuccessRate;

% Information used for confusion matrix calculation and
% the second way to calc mean ROC curve
% meanPScore = meanPScore + PScore;

lenPScore = length(PScore);
for c = 1:nClass
    cl = num2str(c);
    % Receiver operating characteristic (ROC)
    [rocX,rocY,rocT,AUC,OPTROCPT] = perfcurve(targetClassS,...
        PScore(:,c),c,'NBoot',1000,'TVals',[0:1/(lenPScore-1):1],...
        'Options', parOptions);
    % TODO: Save confidence interval of ROC, not only mean value.
    rocData{c,1} = rocX(:,1);
    rocData{c,2} = rocY(:,1);
    rocData{c,3} = rocT;
    rocData{c,4} = AUC(:,1);
    rocData{c,5} = OPTROCPT;
    
    % Class filter (select actual class)
    ii = find(ismember(targetClassS,cl));
    
    % Store class label stats compared with actual class
    for c1=1:nClass
        aux = strcmp(label(ii),int2str(c1));
        classificationStats(c,c1) = classificationStats(c,c1) + sum(aux);
    end
    
    % Overall class true positive (%)
    aux = strcmp(targetClassS(ii),label(ii));
    TP = sum(aux);
    statisticalM(c,1) = 100*TP/length(aux) ...
        + statisticalM(c,1);
    
    % True positive
    statisticalM(c,2) = statisticalM(c,2) + TP;
    
    % False positive
    ii = find(ismember(label,cl));
    aux = ~strcmp(targetClassS(ii),label(ii));
    FP = sum(aux);
    statisticalM(c,3) = statisticalM(c,3) + FP;
    
    % True negative
    ii = find(~ismember(targetClassS,cl));
    aux = ~strcmp(label(ii),cl);
    TN = sum(aux);
    statisticalM(c,4) = statisticalM(c,4) + TN;
    
    % False negative
    ii = find(ismember(targetClassS,cl));
    aux = ~strcmp(label(ii),cl);
    FN = sum(aux);
    statisticalM(c,5) = statisticalM(c,5) + FN;
    
    % Sensitivity
    sensitivity = TP / (TP + FN);
    statisticalM(c,6) = statisticalM(c,6) + sensitivity;
    
    % Specificity
    specificity = TN / (TN + FP);
    statisticalM(c,7) = statisticalM(c,7) + specificity;
    
    % Precision
    precision = TP / (TP + FP);
    statisticalM(c,8) = statisticalM(c,8) + precision;
    
    % Accuracy
    accuracy = (TP + TN) / (TP + TN + FP + FN);
    statisticalM(c,9) = statisticalM(c,9) + accuracy;
    
end
fprintf('------\nStats calcs finished.\n');
%% Print classification statistics (cross-validation)
fprintf(baseName);
fprintf('\nOverall success rate: %.3f%%', overallSuccessRate);
% Save overallSuccessRate in ExperimentsConfig xls file
if useConfigData == 1
    cvmAcuracyCol = colHeadingsLetters{ismember(colHeadingsTxt,'CVMAcc')};
    xlwrite(configFilePath, overallSuccessRate, 1, [cvmAcuracyCol int2str(mainCounter+1)]);
end
% Classification loss for observations not used for training
kFLoss = kfoldLoss(CVMdl, 'Options', parOptions);
fprintf(...
    '\nClassification loss for observations not used for training: %.3f \n\n',...
    kFLoss);
for c = 1:nClass
    fprintf('Class %d overall success rate: %.3f%% \n', c, statisticalM(c,1));
    fprintf('Class %d mean sensitivity: %.3f \n', c, statisticalM(c,6));
    fprintf('Class %d mean specificity: %.3f \n', c, statisticalM(c,7));
    fprintf('Class %d mean precision: %.3f \n', c, statisticalM(c,8));
    fprintf('Class %d mean accuracy: %.3f \n', c, statisticalM(c,9));
    fprintf('Class %d mean ROC AUC: %.3f \n', c, rocData{c,4});
    fprintf('\n');
end

% Close parallel environment (or keep the pool open to save time)
% By default, as setup in Matlab, the parallel pool will be deleted after
% 30 minutes if it is idle.
% delete pool;

%% Generate table to save
measures = {'OverallSuccessRate',...
    'TruePositiveMean','FalsePositiveMean','TrueNegativeMean','FalseNegativeMean',...
    'MeanSensitivity','MeanSpecificity','MeanPrecision','MeanAccuracy'};
% Classification statistical table
CST = array2table(statisticalM,...
    'RowNames', classNames,...
    'VariableNames', measures);
% ROC
rocMeasures = {'MeanX','MeanY','MeanThreshold','MeanAUC','MeanOPTROCPT'};
ROC = cell2struct(rocData, rocMeasures, 2);
if saveData == 'y'
    statFileName = [baseName '_CrossValidationStats.mat'];
    save([iterationPath '/' statFileName],'overallSuccessRate','kFLoss','CST','ROC');
end
%% ROC curve plot and Optimal operating point of the ROC curve
fixROCToPlot = 1;
baseFontSize = 16;
markerSize = 10;
markerType = {'^','*','s','+'};
lineStyle = {'-',':','-.','--'};
% Line style and color
% lineSC = {['k' lineStyle{1} markerType{1}],...
%     ['k' lineStyle{2} markerType{2}],...
%     ['k' lineStyle{3} markerType{3}],...
%     ['k' lineStyle{4} markerType{4}]};
lineSC = {['k' lineStyle{1}],...
    ['k' lineStyle{2}],...
    ['k' lineStyle{3}],...
    ['k' lineStyle{4}]};
figure;
hold on;
iic = cell(nClass,2);
for c = 1:nClass
    if fixROCToPlot == 1
        rocData{c,1}(1) = 0;
        rocData{c,2}(1) = 0;
    end
    plot(rocData{c,1}, rocData{c,2}, lineSC{c}, 'MarkerSize', markerSize);
    aucTaskNames{c} = [groupsNames{c} ' (AUC = ' num2str(rocData{c,4},2) ')'];
    %     Used by follow for, to plot points under original ROC curve
    %     endiic = find(rocData{c,1} == max(rocData{c,1}),1);
    %     iic{c,1} = rocData{c,1}(1:floor(endiic/4):endiic);
    %     iic{c,2} = rocData{c,2}(1:floor(endiic/4):endiic);
end
legend(aucTaskNames, 'Location','Best','fontsize',baseFontSize);
for c = 1:nClass
    %     % Plot points under original ROC curve
    %     plot(iic{c,1}, iic{c,2}, ['k' markerType{c}], 'MarkerSize',markerSize);
    % Plot OPT ROC point
    plot(rocData{c,5}(1), rocData{c,5}(2),['r' markerType{c}],...
        'MarkerSize',markerSize);
end
xlabel('False positive rate','fontsize',baseFontSize);
ylabel('True positive rate','fontsize',baseFontSize);
title(['SVM ROC Curve ' titleBaseName],'fontsize',baseFontSize+2);
plot([0, 1], [0, 1], '--','Color',[0.5,0.5,0.5]);
set(gca,'fontsize',baseFontSize);
ylim([0 1.05]);
xlim([-0.05 1]);
grid on;
hold off;
if saveData == 'y'
    figureName = [baseName '_ROC.fig'];
    savefig([iterationPath '/' figureName]);
    close;
end
%% Confusion Matrix plot
% Plot tips:
% https://www.mathworks.com/matlabcentral/answers/76739-how-can-i-change-color-and-font-size-in-plotconfusion-figures
baseFontSize = 16;
targets = zeros(length(classNumbers),length(x));
x_t = x';
for j = 1:length(targets)
    i = x_t(3,j);
    targets(i,j) = 1;
end
PScore_t_matrix = [];
targetsMatrix = [];
% for i=1:length(pscorePerModel)
%     PScore_t = pscorePerModel{i}';
PScore_t = PScore';
PScore_t_matrix = [PScore_t_matrix PScore_t];
targetsMatrix = [targetsMatrix targets];
% end
figure;
hold on;
pc = plotconfusion(targetsMatrix,PScore_t_matrix);
title(['SVM Confusion Matrix ' titleBaseName]...
    ,'FontSize',baseFontSize);
set(gca,'fontsize',baseFontSize+2);
set(gca,'xticklabel',{groupsNames{:} ''});
set(gca,'yticklabel',{groupsNames{:} ''});
set(findobj(gca,'type','text'),'fontsize',baseFontSize);
hold off;
if saveData == 'y'
    figureName = [baseName '_ConfusionMatrix.fig'];
    savefig([iterationPath '/' figureName]);
    close;
end
%% Save classification statistics in .XLSX (Cross-Validated)

% cv = Confusion value = fraction of samples misclassified
% cm = S-by-S confusion matrix, where cm(i,j) is the number of samples whose target is the ith class that was classified as j
[cv,cm,ind,per] = confusion(targetsMatrix,PScore_t_matrix);
cv = cv * 100;
% Number of samples per class
nSamples = sum(cm,2);
cmP = [];
for i=1:length(nSamples)
    rowCmP = arrayfun(@(x) x/nSamples(i), cm(i,:), ...
        'UniformOutput', false);
    % Percentage normalized between 0 - 1
    cmP = [cmP;rowCmP];
end

% Classification stats (cross-validation) by task and group location at excel file
taskStartRange = {'J6','J12','J18','J24'};
% Classification loss (cross-validation) for data not used for training
lossStartRange = {'N6','N12','N18','N24'};

xlsPath = [workFolder '/Results/'];
if exist('resultsXlsName','var') == 1
    xlsFileName = ['ResultsMar2017-' computeMapAlgLabels{computeMapAlg} ...
        resultsXlsName '.xlsx'];
else
    xlsFileName = ['ResultsMar2017-' computeMapAlgLabels{computeMapAlg} ...
        '.xlsx'];
end

% Check if results xls file already exists, if not, creates one using the
% empty xls file
if size(dir([xlsPath xlsFileName]),1) == 0 && saveData == 'y'
    emptyResultsFilePath = [xlsPath 'ResultsMar2017-Empty.xlsx'];
    copyfile(emptyResultsFilePath, [xlsPath xlsFileName]);
end

sheetName = selectedDataSet;
if saveData == 'y'
    xlwrite([xlsPath xlsFileName], cmP, sheetName, taskStartRange{task});
    xlwrite([xlsPath xlsFileName], cv, sheetName, lossStartRange{task});
    % Save ResultsXLSFileName in ExperimentsConfig xls file
    if useConfigData == 1
        resultsXLSFN = colHeadingsLetters{ismember(colHeadingsTxt,'ResultsXLSFN')};
        xlwrite(configFilePath, {xlsFileName}, 1, [resultsXLSFN int2str(mainCounter+1)]);
    end
end

%% Test set classification

% Train SVM model to predict new data (test set)
% Fix this to use the same model used in cross-validation step
Mdl = fitcecoc(x(:,1:2),targetClassS, 'Learners', template, ...
    'FitPosterior',1, ...
    'ClassNames', classNames, 'Verbose', 1,...
    'Coding', codingDesign, 'Options', parOptions);

% if computeMapAlg ~= 1
%     % Reduce dimensionality of new data using trained model (that learned
%     % Sammon's or t-SNE function)
%     oos.projection = trainedN(bestNN).net(testDataPCA')';
% elseif computeMapAlg == 1
%     oos.projection = testDataPCA;
% end
% Predict new reduced data using trained svm model
[testLabelSVM,~,testPScore,testPosterior] = ...
    predict(Mdl, oos.projection, 'Options', parOptions);
% Classification loss
tsLoss = loss(Mdl, oos.projection, testLabels, 'Options', parOptions);

testTargetLabels = num2cell(testLabels);

% Organize results to be used by confusion functions
testTargets = zeros(length(classNumbers),length(testTargetLabels));
testLabels_t = testLabels';
for j = 1:length(testTargets)
    i = testLabels_t(1,j);
    testTargets(i,j) = 1;
end
testPScore_t_matrix = [];
testTargetsMatrix = [];
testPScore_t = testPScore';

testPScore_t_matrix = [testPScore_t_matrix testPScore_t];
testTargetsMatrix = [testTargetsMatrix testTargets];

% TODO: Compute ROC for test set

if saveData == 'y'
    tsModelFileName = [baseName '_SvmModel.mat'];
    save([iterationPath '/' tsModelFileName],'Mdl','tsLoss',...
        'testPScore_t_matrix','testTargetsMatrix');
end

%% Print classification statistics (test set)
fprintf('\n');
tsOverralSucessRate = 100 - (tsLoss*100);
fprintf('\nOverall test set success rate: %.3f%%', tsOverralSucessRate);
% Save tsOverallSuccessRate in ExperimentsConfig xls file
if useConfigData == 1
    testSetAccCol = colHeadingsLetters{ismember(colHeadingsTxt,'TestSetAcc')};
    xlwrite(configFilePath, tsOverralSucessRate, 1, [testSetAccCol int2str(mainCounter+1)]);
end
fprintf(...
    '\nTest set classification loss: %.3f \n\n',...
    tsLoss);
fprintf('\n');

%% Test set classification statistics (to be saved in .XLSX)

% cv = Confusion value = fraction of samples misclassified
% cm = S-by-S confusion matrix, where cm(i,j) is the number of samples
% whose target is the ith class that was classified as j
[tsCv,tsCm,tsInd,tsPer] = confusion(testTargetsMatrix,testPScore_t_matrix);
tsCv = tsCv * 100;
% Number of samples per class
tsNSamples = sum(tsCm,2);
tsCmP = [];
for i=1:length(tsNSamples)
    tsRowCmP = arrayfun(@(x) x/tsNSamples(i), tsCm(i,:), ...
        'UniformOutput', false);
    % Percentage normalized between 0 - 1
    tsCmP = [tsCmP;tsRowCmP];
end

% Classification stats (test set) by task and group location at excel file
tsTaskStartRange = {'C6','C12','C18','C24'};
% Classification loss (test set) for data not used for training
tsLossStartRange = {'G6','G12','G18','G24'};

if saveData == 'y'
    xlwrite([xlsPath xlsFileName], tsCmP, sheetName, tsTaskStartRange{task});
    xlwrite([xlsPath xlsFileName], tsCv, sheetName, tsLossStartRange{task});
end

%% Confusion matrix plot (test set)
% figure;
% hold on;
% pc = plotconfusion(testTargetsMatrix,testPScore_t_matrix);
% title(['SVM Confusion Matrix (test set) ' titleBaseName]...
%     ,'FontSize',baseFontSize);
% set(gca,'fontsize',baseFontSize+2);
% set(gca,'xticklabel',{groupsNames{:} ''});
% set(gca,'yticklabel',{groupsNames{:} ''});
% set(findobj(gca,'type','text'),'fontsize',baseFontSize);
% hold off;












