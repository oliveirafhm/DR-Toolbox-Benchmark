% By F?bio Henrique (oliveirafhm@gmail.com)
% Run at the end of everything
% ...
% 06/2017

%% Load results from ResultsControl xlsx file

[resultsControlFileName, resultsControlPathName] = uigetfile('.xlsx', ...
    'Select results control xlsx file');
resultsControlFilePath = strcat(resultsControlPathName, resultsControlFileName);

[~,~,resultsControl] = xlsread(resultsControlFilePath,1);
resultsControlHeadings = resultsControl(1,1:end);
resultsControlData = resultsControl(2:end,1:end);

%% Organize ResultsControl data
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
yCV = NaN(length(resultsControlData), length(DRAlgs), nGroups);
yT = NaN(length(resultsControlData), length(DRAlgs), nGroups);
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
plotBar = 0;

nBars = nGroups * nSets; % length(GMCols)/2
meanGM = zeros(length(DRAlgs),nBars);
meanSTD = zeros(length(DRAlgs),nBars);
%
clearvars GM GM_CV GM_T;

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
if plotBar
    %barGroups = {'S_{H}','S_{PD}','S_{DBS}'};
    legend = {'PCA', 'Sammon''s mapping', ['t-SNE ' tSneL]};
    
    % Mean GM chart of training set with LOO cross-validation
    barY1 = meanGM(:,1:3)';
    barY1STD = meanSTD(:,1:3)';
    figure;
    barweb(barY1, barY1STD, [], groups, ...
        'True positive general mean for each group of subjects of training set (LOO CV)',...
        'Groups of subjects', 'Success rate (%)', bone, 'y', legend, 2, 'plot', 10);
    
    % Mean GM chart of test set
    barY2 = meanGM(:,4:6)';
    barY2STD = meanSTD(:,4:6)';
    figure;
    barweb(barY2, barY2STD, [], groups, ...
        'True positive general mean for each group of subjects of test set',...
        'Groups of subjects', 'Success rate (%)', bone, 'y', legend, 2, 'plot', 10);
end
%% Boxplot for each DR method to show each group
% mean (Grand boxplot). For LOO CV and Test set.
legend = {'PCA', 'Sammon''s mapping', ['t-SNE ' tSneL]};
labels = {'PCA', 'Sammon''s', 't-SNE '};
factor = 10;

%https://github.com/IoSR-Surrey/MatlabToolbox
boxplot2 = @iosr.statistics.boxPlot;

% LOO CV TP boxplot
figure;
boxplot2(groups, yCV, 'notch', true, ...
    'symbolColor','k',...
    'medianColor','k',...
    'style','hierarchy',...
    'xSeparator',true,...
    'groupLabels',{labels},...
    'groupLabelFontSize', 9+factor,...
    'theme', 'colorboxes',...
    'outlierSize', 36+factor);
% ylim([round(min(min(min(yCV))),1) - .05 round(max(max(max(yCV))),1)]);
box on;
set(gca,'ygrid','on', 'FontSize',11+factor*0.9);
title('True positive general mean distribution for each group of subjects (cross-validation step)',...
    'FontSize',11+factor*1.3);
set(findobj(gca,'Type','text'),'FontSize',11+factor*0.9);
xlabel('','FontSize', 11+factor);
ylabel('Success rate (%)', 'FontSize', 11+factor);

% Test TP boxplot
figure;
boxplot2(groups, yT, 'notch', true, ...
    'symbolColor','k',...
    'medianColor','k',...
    'style','hierarchy',...
    'xSeparator',true,...
    'groupLabels',{labels},...
    'groupLabelFontSize', 9+factor,...
    'theme', 'colorboxes',...
    'outlierSize', 36+factor);
ylim([round(min(min(min(yT))),1) - .05 round(max(max(max(yT))),1)]);
box on;
set(gca,'ygrid','on', 'FontSize',11+factor*0.9);
title('True positive general mean distribution for each group of subjects (test step)',...
    'FontSize',11+factor*1.3);
set(findobj(gca,'Type','text'),'FontSize',11+factor*0.9);
xlabel('','FontSize', 11+factor);
ylabel('Success rate (%)', 'FontSize', 11+factor);

%% Overall confusion matrix of general mean (Grand average)
workFolder = 'PaperNeural';
resultsPath = [workFolder '/Results/'];
saveGrandConfusion = 'y';
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
             oAccCV / sum(sum(grandConfusionMT(1:nGroups,1:nGroups,drA)))*100;
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
