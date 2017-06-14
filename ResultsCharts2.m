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
    for c = 1:length(DRAlgs)
        if c ~= 3
            idx{c} = find([resultsControlData{:,2}] == c);
        else
            idx0 = find([resultsControlData{:,2}] == c);
            % Split XLSX file name
            subStrings = cellfun(@(s) strsplit(s,'_'), ...
                resultsControlData(idx0,xlsFileNameCol), ...
                'UniformOutput', false);  
            idx{c} = [];
            for ssi = 1:length(subStrings)
                % Get perplexity value from Xls file name string
                ss = strsplit(subStrings{ssi}{end},'.');
                ss = str2double(ss{1});
                if ss > perpThreshold(1) && ss < perpThreshold(2)
                    idx{c}(end+1) = idx0(ssi);
                end
            end
            tSneL = ['(' int2str(perpThreshold(1)) ' < perp < ' ...
                int2str(perpThreshold(2)) ')'];
        end
        % row = method accuracy, column = group of subject, page = DR method
        for r = 1:length(idx{c})
            %             yCV(r,c,p) = [resultsControlData{idx(r),GM_CV_Cols(p)}];
            %             yT(r,c,p) = [resultsControlData{idx(r),GM_T_Cols(p)}];
            yCV(r,p,c) = [resultsControlData{idx{c}(r),GM_CV_Cols(p)}];
            yT(r,p,c) = [resultsControlData{idx{c}(r),GM_T_Cols(p)}];
        end
    end
end
%% Bar chart with std (error bars) for each DR method to show each group
% mean and std. For LOO CV and Test set.
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
% mean. For LOO CV and Test set.
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

%% Overall confusion matrix of general mean


