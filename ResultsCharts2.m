% By F?bio Henrique (oliveirafhm@gmail.com)
% Run at the end of everything
% ...
% 06/2017 (last modification: 04/01/2018)
%% ------------------------ Part 1 ------------------------
poiPath = 'Third party codes/xlswrite/20130227_xlwrite/poi_library/';
javaaddpath([poiPath 'poi-3.8-20120326.jar']);
javaaddpath([poiPath 'poi-ooxml-3.8-20120326.jar']);
javaaddpath([poiPath 'poi-ooxml-schemas-3.8-20120326.jar']);
javaaddpath([poiPath 'xmlbeans-2.3.0.jar']);
javaaddpath([poiPath 'dom4j-1.6.1.jar']);
javaaddpath([poiPath 'stax-api-1.0.1.jar']);
addpath(genpath('Third party codes'));

%% Load results from ResultsControl xlsx file
% ------------------------ Part 1 ------------------------

[resultsControlFileName, resultsControlPathName] = uigetfile('.xlsx', ...
    'Select results control xlsx file');
resultsControlFilePath = strcat(resultsControlPathName, resultsControlFileName);

[~,~,resultsControl] = xlsread(resultsControlFilePath,1);
resultsControlHeadings = resultsControl(1,1:end);
resultsControlData = resultsControl(2:end,1:end);

%% Organize ResultsControl data
% ------------------------ Part 1 ------------------------
% Columns with general mean and std
GMCols = 4:2:length(resultsControlHeadings);
GM_CV_Cols = GMCols(1:length(GMCols)/2);
GM_T_Cols = GMCols(length(GMCols)/2+1:end);
GMSTDCols = 5:2:length(resultsControlHeadings);
GMSTD_CV_Cols = GMSTDCols(1:length(GMSTDCols)/2);
GMSTD_T_Cols = GMSTDCols(length(GMSTDCols)/2+1:end);
% N of DR algorithms
DRAlgs = unique(cell2mat(resultsControlData(:,2)));
% Sh , Spd, Sdbs
nGroups = 3;
groups = {'S_{H}','S_{PD}','S_{DBS}'};
% N of sets: LOO CV set , test set
nSets = 2;
%Only for t-SNE
perpThreshold = [17,51];
xlsFileNameCol = 3;
% ResultsControl data organized in an 3D matrix
yCV = NaN(length(resultsControlData), nGroups, length(DRAlgs));
yT = NaN(length(resultsControlData), nGroups, length(DRAlgs));
idx = {};
for p = 1:nGroups
    for j = 1:length(DRAlgs)
        if j ~= 3
            idx{j} = find([resultsControlData{:,2}] == j);
        else
            idx0 = find([resultsControlData{:,2}] == j);
            % Split XLSX file name
            subStrings = cellfun(@(s) strsplit(s,'_'), ...
                resultsControlData(idx0,xlsFileNameCol), ...
                'UniformOutput', false);
            idx{j} = [];
            for ssi = 1:length(subStrings)
                % Get perplexity value from Xls file name string
                ss = strsplit(subStrings{ssi}{end},'.');
                ss = str2double(ss{1});
                if ss > perpThreshold(1) && ss < perpThreshold(2)
                    idx{j}(end+1) = idx0(ssi);
                end
            end
            tSneL = ['(' int2str(perpThreshold(1)) ' < perp < ' ...
                int2str(perpThreshold(2)) ')'];
        end
        % row = method accuracy, column = group of subject, page = DR method
        for r = 1:length(idx{j})
            %             yCV(r,c,p) = [resultsControlData{idx(r),GM_CV_Cols(p)}];
            %             yT(r,c,p) = [resultsControlData{idx(r),GM_T_Cols(p)}];
            yCV(r,p,j) = [resultsControlData{idx{j}(r),GM_CV_Cols(p)}];
            yT(r,p,j) = [resultsControlData{idx{j}(r),GM_T_Cols(p)}];
        end
    end
end
%% Bar chart with std (error bars) for each DR method to show each group
% mean (grand average) and std. For LOO CV and Test set.
% ------------------------ Part 1 ------------------------
% plotBar = 0;
% 
% nBars = nGroups * nSets; % length(GMCols)/2
% meanGM = zeros(length(DRAlgs),nBars);
% meanSTD = zeros(length(DRAlgs),nBars);
% %
% clearvars GM GM_CV GM_T;

% for j = 1:length(DRAlgs)
% %     idx = find([resultsControlData{:,2}] == j);
%     if j ~= 3
%         idx = find([resultsControlData{:,2}] == j);
%     % Get perplexity value from str to filter
%     else
%         idx0 = find([resultsControlData{:,2}] == j);
%         % Split XLSX file name
%         subStrings = cellfun(@(s) strsplit(s,'_'), ...
%             resultsControlData(idx0,3), 'UniformOutput', false);
%         idx = [];
%         for ssi = 1:length(subStrings)
%             % Get perplexity value from Xls file name string
%             ss = strsplit(subStrings{ssi}{end},'.');
%             ss = str2double(ss{1});
%             if ss > perpThreshold(1) && ss < perpThreshold(2)
%                 idx(end+1) = idx0(ssi);
%             end
%         end
%         tSneL = ['(' int2str(perpThreshold(1)) ' < perp < ' int2str(perpThreshold(2)) ')'];
%     end
%     kk = 1; % used to avoid std col
%     kj = j; % used to sort data for boxplot
%     for k = 1:nBars
%         % Used in bar plot
%         meanGM(j,k) = mean([resultsControlData{idx,GMCols(kk)}]);
%         meanSTD(j,k) = mean([resultsControlData{idx,GMCols(kk+1)}]);
%         % Used in box plot
%         if k <= nBars / 2
%             if exist('GM_CV','var')
%                 GM_CV(kj) = {[resultsControlData{idx,GMCols(kk)}]};
%             else
%                 GM_CV = {[resultsControlData{idx,GMCols(kk)}]};
%             end
%         else
%             if exist('GM_T','var')
%                 GM_T(kj) = {[resultsControlData{idx,GMCols(kk)}]};
%             else
%                 GM_T = {[resultsControlData{idx,GMCols(kk)}]};
%             end
%         end
%         kk = kk + 2;
%         kj = kj + length(DRAlgs);
%         % Reset kj for GM_T sorting
%         if k == nBars/2, kj = j; end
%     end
% end
% if plotBar
%     %barGroups = {'S_{H}','S_{PD}','S_{DBS}'};
%     legend = {'PCA', 'Sammon''s mapping', ['t-SNE ' tSneL]};
%     
%     % Mean GM chart of training set with LOO cross-validation
%     barY1 = meanGM(:,1:3)';
%     barY1STD = meanSTD(:,1:3)';
%     figure;
%     barweb(barY1, barY1STD, [], groups, ...
%         'True positive general mean for each group of subjects of training set (LOO CV)',...
%         'Groups of subjects', 'Success rate (%)', bone, 'y', legend, 2, 'plot', 10);
%     
%     % Mean GM chart of test set
%     barY2 = meanGM(:,4:6)';
%     barY2STD = meanSTD(:,4:6)';
%     figure;
%     barweb(barY2, barY2STD, [], groups, ...
%         'True positive general mean for each group of subjects of test set',...
%         'Groups of subjects', 'Success rate (%)', bone, 'y', legend, 2, 'plot', 10);
% end
%% Calc mean of grand average success rate
nanmean(yCV)
nanmean(yT)

%% Boxplot for each DR method to show each group
% mean (Grand boxplot). For LOO CV and Test set.
% ------------------------ Part 1 ------------------------
legend = {'PCA', 'Sammon''s mapping', ['t-SNE ' tSneL]};
labels = {'PCA', 'Sammon''s', 't-SNE '};
factor = 10;

%https://github.com/IoSR-Surrey/MatlabToolbox
boxplot2 = @iosr.statistics.boxPlot;

% LOO CV TP boxplot
figure;
boxplot2(groups, yCV, 'notch', false, ...
    'symbolColor','k',...
    'medianColor','k',...
    'style','hierarchy',...
    'xSeparator',true,...
    'groupLabels',{labels},...
    'groupLabelFontSize', 9+factor,...
    'theme', 'default',...
    'outlierSize', 36+factor);
% ylim([round(min(min(min(yCV))),1) - .05 round(max(max(max(yCV))),1)]);
box on;
set(gca,'ygrid','on', 'FontSize',11+factor*0.9);
title('True positive general mean distribution for each group of subjects (cross-validation step)',...
    'FontSize',11+factor*1.3);
set(findobj(gca,'Type','text'),'FontSize',11+factor*0.9);
xlabel('','FontSize', 11+factor);
ylabel('Success rate (%)', 'FontSize', 11+factor);
pbaspect([1.5 1 1]); % Figure aspect ratio

% Test TP boxplot
figure;
boxplot2(groups, yT, 'notch', false, ...
    'symbolColor','k',...
    'medianColor','k',...
    'style','hierarchy',...
    'xSeparator',true,...
    'groupLabels',{labels},...
    'groupLabelFontSize', 9+factor,...
    'theme', 'default',...
    'outlierSize', 36+factor);
ylim([round(min(min(min(yT))),1) - .05 round(max(max(max(yT))),1)]);
box on;
set(gca,'ygrid','on', 'FontSize',11+factor*0.9);
title('True positive general mean distribution for each group of subjects (test step)',...
    'FontSize',11+factor*1.3);
set(findobj(gca,'Type','text'),'FontSize',11+factor*0.9);
xlabel('','FontSize', 11+factor);
ylabel('Success rate (%)', 'FontSize', 11+factor);
pbaspect([1.5 1 1]); % Figure aspect ratio
%% Normality test
% [h,p] = kstest(yCV(:,3,3));

%% Statistical analysis of groups in boxplot
% ------------------------ Part 1 ------------------------
savePValues = 'n';
% DRAlgs X DRAlgs X Group of subjects
pValueCV = NaN(length(DRAlgs), length(DRAlgs), nGroups);
pValueT = NaN(length(DRAlgs), length(DRAlgs), nGroups);

for g = 1:nGroups
    for pi = 1:length(DRAlgs)
        for pj = 1:length(DRAlgs)
            %             [pCV,hCV] = ranksum(yCV(:,g,pi),yCV(:,g,pj));
            [hCV,pCV] = kstest2(yCV(:,g,pi),yCV(:,g,pj));
            pValueCV(pi, pj, g) = pCV;
            
            %             [pT,hT] = ranksum(yT(:,g,pi),yT(:,g,pj));
            [hT,pT] = kstest2(yT(:,g,pi),yT(:,g,pj));
            pValueT(pi, pj, g) = pT;
        end
    end
end

% Save p-value stats of true positive general mean in ResultsControlFinal.xlsx file
if savePValues == 'y'
    sheetNumber = 4;% hard code
    testStartRange = {'C7','C15','C23'};
    crossValidatedStartRange = {'J7','J15','J23'};
    for i=1:nGroups
        xlwrite(resultsControlFilePath, pValueCV(:,:,i),...
            sheetNumber, crossValidatedStartRange{i});
        xlwrite(resultsControlFilePath, pValueT(:,:,i),...
            sheetNumber, testStartRange{i});
    end
end
%% Overall confusion matrix of general mean (Grand average)
% ------------------------ Part 1.1 ------------------------
workFolder = 'PaperNeural';
resultsPath = [workFolder '/Results/'];
saveGrandConfusion = 'n';
grandConfusionMCV = zeros(nGroups+1, nGroups+1, length(DRAlgs));
grandConfusionMT = zeros(nGroups+1, nGroups+1, length(DRAlgs));
for drA = 1:length(DRAlgs)
    for i=1:length(idx{drA})
        % Load group experiment xls file (ResultsMar2017..) - Summary sheet
        % For Cross-Validation
        [~,~,resultsSumRange] = xlsread([resultsPath...
            resultsControlData{idx{drA}(i),xlsFileNameCol}], 1, 'Q73:Y73');
        resultsSumRange = [resultsSumRange{:}];
        % Construct grand confusion 3d matrix for CV
        gc = 1;
        for j=1:nGroups:nGroups*nGroups
            gl = 1;
            for jj=j:j+nGroups-1
                grandConfusionMCV(gl, gc, resultsControlData{idx{drA}(i),2}) = ...
                    grandConfusionMCV(gl, gc, resultsControlData{idx{drA}(i),2}) + ...
                    resultsSumRange(jj);
                gl = gl + 1;
            end
            gc = gc + 1;
        end
        % Load group experiment xls file (ResultsMar2017..) - Summary sheet
        % For Test set
        [~,~,resultsSumRange] = xlsread([resultsPath...
            resultsControlData{idx{drA}(i),xlsFileNameCol}], 1, 'C73:K73');
        resultsSumRange = [resultsSumRange{:}];
        % Construct grand confusion 3d matrix for test set
        gc = 1;
        for j=1:nGroups:nGroups*nGroups
            gl = 1;
            for jj=j:j+nGroups-1
                grandConfusionMT(gl, gc, resultsControlData{idx{drA}(i),2}) = ...
                    grandConfusionMT(gl, gc, resultsControlData{idx{drA}(i),2}) + ...
                    resultsSumRange(jj);
                gl = gl + 1;
            end
            gc = gc + 1;
        end
    end
    % Calc average of grand confusion 3d matrix
    grandConfusionMCV(:,:,drA) = grandConfusionMCV(:,:,drA)./length(idx{drA});
    grandConfusionMT(:,:,drA) = grandConfusionMT(:,:,drA)./length(idx{drA});
    % Calc last column, last row and last cell of confusion matrix
    oAccCV = 0; oAccT = 0;
    for gIdx = 1:nGroups
        grandConfusionMCV(gIdx,nGroups+1,drA) = ...
            grandConfusionMCV(gIdx,gIdx,drA)*100/sum(grandConfusionMCV(gIdx,:,drA));
        grandConfusionMCV(nGroups+1,gIdx,drA) = ...
            grandConfusionMCV(gIdx,gIdx,drA)*100;
        grandConfusionMT(gIdx,nGroups+1,drA) = ...
            grandConfusionMT(gIdx,gIdx,drA)*100/sum(grandConfusionMT(gIdx,:,drA));
        grandConfusionMT(nGroups+1,gIdx,drA) = ...
            grandConfusionMT(gIdx,gIdx,drA)*100;
        % Calc of overall accuracy
        oAccCV = oAccCV + grandConfusionMCV(gIdx,gIdx,drA);
        oAccT = oAccT + grandConfusionMT(gIdx,gIdx,drA);
        if gIdx == nGroups
            grandConfusionMCV(end,end,drA) = ...
                oAccCV / sum(sum(grandConfusionMCV(1:nGroups,1:nGroups,drA)))*100;
            grandConfusionMT(end,end,drA) = ...
                oAccT / sum(sum(grandConfusionMT(1:nGroups,1:nGroups,drA)))*100;
        end
    end
end
% Save grand confusion matrix in ResultsControlFinal.xlsx file
if saveGrandConfusion == 'y'
    sheetNumber = 3;% hard code
    testStartRange = {'C6','C14','C22'};
    crossValidatedStartRange = {'J6','J14','J22'};
    for i=1:length(DRAlgs)
        xlwrite(resultsControlFilePath, grandConfusionMT(:,:,i),...
            sheetNumber, testStartRange{i});
        xlwrite(resultsControlFilePath, grandConfusionMCV(:,:,i),...
            sheetNumber, crossValidatedStartRange{i});
    end
end

%% Average ROC with confidence interval for each DR method
% ------------------------ Part 2 ------------------------
% Load experiment config file
[configFileName, configPathName] = uigetfile('.xlsx', ...
    'Select experiment config file');
configFilePath = strcat(configPathName, configFileName);
[~,~,rawConfig] = xlsread(configFilePath,1);
[~,~,colHeadings] = xlsread(configFilePath,2);
colHeadingsTxt = colHeadings(1,:);
colHeadingsLetters = colHeadings(2,:);

expColHeadings = rawConfig(1,:);
expConfigData = cell2struct(rawConfig(2:end,1:end),...
    expColHeadings(1:end), 2);

% Filter DR data
idx = {};
%Only for t-SNE
perpThreshold = [17,51];
nDR = length(unique([expConfigData(:).DRAlg]));
for i = 1:nDR
    if i ~= 3
        idx{i} = find([expConfigData(:).DRAlg] == i);
    else
        % Filter perp of t-sne (like did before for boxplot and confusion matrix)
        idx{i} = find([expConfigData(:).Perplexity] > perpThreshold(1)...
            & [expConfigData(:).Perplexity] < perpThreshold(2));
    end
end

%% ROC calculus for LOOCV
% ------------------------ Part 2 ------------------------
addpath(genpath('PaperNeural'));
parOptions = statset('UseParallel',1);
groupsNames = {'S_{H}','S_{PD}','S_{DBS}'};
% max n of roc samples X n of DR methods x n of groups
rocData = cell(length(idx{3}),nDR,length(groupsNames));
tic
for i = 1:nDR
    for j = 1:length(idx{i})
        % Load: ..._LooSvmModel.mat
        load([expConfigData(idx{i}(j)).BaseFileName '_LooSvmModel.mat']);
        [label,~,PScore,Posterior] = kfoldPredict(CVMdl, 'Options', parOptions);
        lenPScore = length(PScore);
        for c = 1:length(groupsNames)
            [rocX,rocY,rocT,AUC,OPTROCPT] = perfcurve(CVMdl.Y,...
                PScore(:,c),c,'NBoot',1000,'TVals',[0:1/(lenPScore-1):1],...
                'Options', parOptions);
            rocData{j,i,c} = {rocX,rocY,rocT,AUC,OPTROCPT};
        end
        %clear vars
    end
end
toc
% TODO: Save vars
%% ROC calculus for test set (hold-out)
% ------------------------ Part 2 ------------------------
addpath(genpath('PaperNeural'));
parOptions = statset('UseParallel',1);
groupsNames = {'S_{H}','S_{PD}','S_{DBS}'};
% max n of roc samples X n DR methods x n of groups
rocData = cell(length(idx{3}),nDR,length(groupsNames));

tic
for i = 1:nDR
    for j = 1:length(idx{i})
        % Load: ..._OosModelData.mat
        load([expConfigData(idx{i}(j)).BaseFileName '_OosModelData.mat']);
        % Load: ..._SvmModel.mat
        load([expConfigData(idx{i}(j)).BaseFileName '_SvmModel.mat']);
        [testLabelSVM,~,testPScore,testPosterior] = ...
            predict(Mdl, oos.projection, 'Options', parOptions);
        lenPScore = length(testPScore);
        for c = 1:length(groupsNames)
            %             [rocX,rocY,rocT,AUC,OPTROCPT] = perfcurve(testLabels,...
            %                 testPScore(:,c),c,'NBoot',1000,'TVals',[0:1/(lenPScore-1):1],...
            %                 'Options', parOptions);
            [rocX,rocY,rocT,AUC,OPTROCPT] = perfcurve(testLabels,...
                testPScore(:,c),c,'Options', parOptions);
            rocData{j,i,c} = {rocX,rocY,rocT,AUC,OPTROCPT};
        end
        %clear vars
    end
end
toc
% TODO: Save vars
%% Calc mean ROC for LOO
% ------------------------ Part 3 ------------------------
% TODO: Load LOO ROC Data first

% DR algs X Subject group
meanRocDataLoo = cell(nDR, length(groupsNames));
% rocX = []; rocY = []; rocT = []; AUC = []; OPTROCPT = [];
nRocVars = 5;

for c = 1:nDR
    for p = 1:length(groupsNames)
        meanRocDataLoo{c,p} = cell(1,nRocVars);
        for r = 1:length(idx{c})
            for nCell = 1:nRocVars
                if r ~= 1
                    meanRocDataLoo{c,p}{nCell} = ...
                        meanRocDataLoo{c,p}{nCell} + rocData{r,c,p}{nCell};
                else
                    meanRocDataLoo{c,p}{nCell} = rocData{r,c,p}{nCell};
                end
            end
        end
        % Calc mean for each ROC var
        for nCell = 1:nRocVars
            meanRocDataLoo{c,p}{nCell} = ...
                meanRocDataLoo{c,p}{nCell} ./ length(idx{c});
        end
    end
end

%% Mean ROC curve plot for LOOCV
% ------------------------ Part 3 ------------------------
addpath(genpath('Third party codes'));

nClass = length(groupsNames);
drNames = {'PCA', 'Sammon''s', 't-SNE'};

fixROCToPlot = 1;
baseFontSize = 16;

% markerSize = 10;
% markerType = {'^','*','s','+'};
lineStyle = {'-',':','-.','--'};
lineSC = {['k' lineStyle{1}],...
    ['k' lineStyle{2}],...
    ['k' lineStyle{3}],...
    ['k' lineStyle{4}]};
% aucDR = cell([1 7]);
aucDR = {};
for g = 1:nClass
    drRocData = meanRocDataLoo(:,g);
    figure;
    hold on;
    for c = 1:nDR
        rocData = drRocData{c};
        if fixROCToPlot == 1
            rocData{1}(1,1:3) = 0;
            rocData{2}(1,1:3) = 0;
        end
        %         [l,p] =
        boundedline(rocData{1}(:,1), rocData{2}(:,1), ...
            [rocData{2}(:,1)-rocData{2}(:,2) ...
            rocData{2}(:,3)-rocData{2}(:,1)], lineSC{c}, 'alpha');
        %         outlinebounds(l,p);
        %         aucDR{c+c-1} = '';
        %         aucDR{c+c} = [drNames{c} ' (Mean AUC = ' num2str(rocData{4}(1),2) ')'];
        aucDR{c} = [drNames{c} ' (Mean AUC = ' num2str(rocData{4}(1),2) ')'];
    end
    xlabel('False positive rate','fontsize',baseFontSize);
    ylabel('True positive rate','fontsize',baseFontSize);
    title(['Leave-one-out ' groupsNames{g}],'fontsize',baseFontSize+2);
    plot([1, 0], [0, 1], '--','Color',[0.5,0.5,0.5]);
    %     aucDR{end} = '';
    hl = findobj(gca,'Type','line');
    hp = findobj(gca,'Type','patch');
    legend([hl(4) hl(3) hl(2)], aucDR, 'Location','Best','fontsize',baseFontSize);
    set(gca,'fontsize',baseFontSize);
    % ylim([0 1.05]);
    % xlim([-0.05 1]);
    % grid on;
    hold off;
end
pbaspect([1 1 1]);

%% Calc mean ROC for Holdout
% ------------------------ Part 4 ------------------------
% TODO: Load Test set ROC Data first
ciFlag = 0; % confidence interval flag
% DR algs X Subject group
meanRocDataHoldout = cell(nDR, length(groupsNames));
% rocX = []; rocY = []; rocT = []; AUC = []; OPTROCPT = [];
nRocVars = 5;

for c = 1:nDR
    for p = 1:length(groupsNames)
        meanRocDataHoldout{c,p} = cell(1,nRocVars);
        for r = 1:length(idx{c})
            for nCell = 1:nRocVars
                if ciFlag == 0
                    % Fix ROC length to 20 points to enable average
                    n = 20;
                    if nCell == 1 && length(rocData{r,c,p}{nCell}) < n
                        [x,y,ux,iux] = ROCInterp(rocData{r,c,p}{nCell}, ...
                            rocData{r,c,p}{nCell+1}, n, [0 1], '', 0);
                        rocData{r,c,p}{nCell} = x';
                        rocData{r,c,p}{nCell+1} = y';
                        rocData{r,c,p}{nCell+2} = [1:-(1/(n-1)):0]';
                    end
                end
                %
                if r ~= 1
                    meanRocDataHoldout{c,p}{nCell} = ...
                        meanRocDataHoldout{c,p}{nCell} + rocData{r,c,p}{nCell};
                else
                    meanRocDataHoldout{c,p}{nCell} = rocData{r,c,p}{nCell};
                end
            end
        end
        % Calc mean for each ROC var
        for nCell = 1:nRocVars
            meanRocDataHoldout{c,p}{nCell} = ...
                meanRocDataHoldout{c,p}{nCell} ./ length(idx{c});
        end
    end
end

%% Mean ROC curve plot for Holdout
% ------------------------ Part 4 ------------------------
addpath(genpath('Third party codes'));

nClass = length(groupsNames);
drNames = {'PCA', 'Sammon''s', 't-SNE'};

fixROCToPlot = 1;
baseFontSize = 16;

% markerSize = 10;
% markerType = {'^','*','s','+'};
lineStyle = {'-',':','-.','--'};
lineSC = {['k' lineStyle{1}],...
    ['k' lineStyle{2}],...
    ['k' lineStyle{3}],...
    ['k' lineStyle{4}]};
% aucDR = cell([1 7]);
aucDR = {};
for g = 1:nClass
    drRocData = meanRocDataHoldout(:,g);
    figure;
    hold on;
    for c = 1:nDR
        rocData = drRocData{c};
        if fixROCToPlot
            if ciFlag
                rocData{1}(1,1:3) = 0;
                rocData{2}(1,1:3) = 0;
            else
                rocData{1}(1) = 0;
                rocData{2}(1) = 0;
            end
        end
        if ciFlag
            boundedline(rocData{1}(:,1), rocData{2}(:,1), ...
                [rocData{2}(:,1)-rocData{2}(:,2) ...
                rocData{2}(:,3)-rocData{2}(:,1)], lineSC{c}, 'alpha');
        else
            plot(rocData{1}(:,1), rocData{2}(:,1), lineSC{c});
        end
        
        aucDR{c} = [drNames{c} ' (Mean AUC = ' num2str(rocData{4}(1),2) ')'];
    end
    xlabel('False positive rate','fontsize',baseFontSize);
    ylabel('True positive rate','fontsize',baseFontSize);
    title(['Hold-out ' groupsNames{g}],'fontsize',baseFontSize+2);
    if ciFlag
        plot([1, 0], [0, 1], '--','Color',[0.5,0.5,0.5]);
        hl = findobj(gca,'Type','line');
        hp = findobj(gca,'Type','patch');
        legend([hl(4) hl(3) hl(2)], aucDR, 'Location','Best',...
            'fontsize',baseFontSize);
    else
        legend(aucDR, 'Location','Best','fontsize',baseFontSize);
        plot([1, 0], [0, 1], '--','Color',[0.5,0.5,0.5]);
    end
    set(gca,'fontsize',baseFontSize);
    hold off;
end
pbaspect([1 1 1]);