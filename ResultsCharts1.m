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

%% Bar chart with std (error bars) for each DR method to show each group
% mean and std. For LOO CV and Test set.
plotBar = 0;
% Number of columns with general mean and std
GMCols = 4:length(resultsControlHeadings);
DRAlgs = unique(cell2mat(resultsControlData(:,2)));
% Sh , Spd, Sdbs
nGroups = 3;
% N of sets: LOO CV set , test set
nSets = 2;
nBars = nGroups * nSets; % length(GMCols)/2
meanGM = zeros(length(DRAlgs),nBars);
meanSTD = zeros(length(DRAlgs),nBars);
%Only for t-SNE
perpThreshold = [17,51];
% 
clearvars GM GM_CV GM_T;
for j = 1:length(DRAlgs)
%     idx = find([resultsControlData{:,2}] == j);
    if j ~= 3 
        idx = find([resultsControlData{:,2}] == j);
    % Get perplexity value from str to filter
    else        
        idx0 = find([resultsControlData{:,2}] == j);
        % Split XLSX file name
        subStrings = cellfun(@(s) strsplit(s,'_'), ...
            resultsControlData(idx0,3), 'UniformOutput', false);
        idx = [];
        for ssi = 1:length(subStrings)
            % Get perplexity value from Xls file name string
            ss = strsplit(subStrings{ssi}{end},'.');
            ss = str2double(ss{1});
            if ss > perpThreshold(1) && ss < perpThreshold(2)
                idx(end+1) = idx0(ssi); 
            end
        end
        tSneL = ['(' int2str(perpThreshold(1)) ' < perp < ' int2str(perpThreshold(2)) ')'];
    end
    kk = 1; % used to avoid std col
    kj = j; % used to sort data for boxplot
    for k = 1:nBars
        % Used in bar plot
        meanGM(j,k) = mean([resultsControlData{idx,GMCols(kk)}]);
        meanSTD(j,k) = mean([resultsControlData{idx,GMCols(kk+1)}]);
        % Used in box plot
        if k <= nBars / 2
            if exist('GM_CV','var')
                GM_CV(kj) = {[resultsControlData{idx,GMCols(kk)}]};
            else
                GM_CV = {[resultsControlData{idx,GMCols(kk)}]};
            end
        else
            if exist('GM_T','var')
                GM_T(kj) = {[resultsControlData{idx,GMCols(kk)}]};
            else
                GM_T = {[resultsControlData{idx,GMCols(kk)}]};
            end
        end
        kk = kk + 2;
        kj = kj + length(DRAlgs);
        % Reset kj for GM_T sorting
        if k == nBars/2, kj = j; end
    end
end
if plotBar
    barGroups = {'S_{H}','S_{PD}','S_{DBS}'};
    legend = {'PCA', 'Sammon''s mapping', ['t-SNE ' tSneL]};
    
    % Mean GM chart of training set with LOO cross-validation
    barY1 = meanGM(:,1:3)';
    barY1STD = meanSTD(:,1:3)';
    figure;
    barweb(barY1, barY1STD, [], barGroups, ...
        'True positive general mean for each group of subjects of training set (LOO CV)',...
        'Groups of subjects', 'Success rate (%)', bone, 'y', legend, 2, 'plot', 10);
    
    % Mean GM chart of test set
    barY2 = meanGM(:,4:6)';
    barY2STD = meanSTD(:,4:6)';
    figure;
    barweb(barY2, barY2STD, [], barGroups, ...
        'True positive general mean for each group of subjects of test set',...
        'Groups of subjects', 'Success rate (%)', bone, 'y', legend, 2, 'plot', 10);
end
%% Boxplot for each DR method to show each group
% mean. For LOO CV and Test set.

% Boxplot re-wrapping to enable the use of cell array containing vectors of
% variable size
col=@(x)reshape(x,numel(x),1);
boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),...
    cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});

%boxplot2 = @iosr.statistics.boxPlot;

% boxplot2 automatically generates the necessary grouping array. 
% All you need to pass to it is a cell array of the vectors you want box plotted.
% Use GM var from last section
GM = [GM_CV,GM_T];
figure;
boxplot2(GM', 'Notch','true');
set(gca,'ygrid','on');
title('True positive general mean distribution for each group of subjects');


%% Line chart to compare overall results of each method (PCA, Sammon, t-SNE)
