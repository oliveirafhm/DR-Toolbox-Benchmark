% By F?bio Henrique (oliveirafhm@gmail.com)
% Run after whole experiments
% Fit experiments data in ResultsControl xls file
% 06/2017
%%
% Load ExperimentsConfig xls file again (now with all columns filled)
% Insert alert that this script uses vars from anothers scripts
% configFilePath, ...
[~,~,rawConfig2] = xlsread(configFilePath,1);
expColHeadings2 = rawConfig2(1,1:end);
expConfigData2 = rawConfig2(2:end,1:end);
% Counts how many experiments were finished
colHeadingsN = 1:length(colHeadingsLetters);
resultsXLSFNColNumber = colHeadingsN(ismember(colHeadingsTxt,'ResultsXLSFN'));
%resultsXLSFNColNumber = 21; % or last one or T column
% If the experiment was not completed, we replace NaN by empty string 
for i = 1:length(expConfigData2(:,resultsXLSFNColNumber))
    if isnan(expConfigData2{i,resultsXLSFNColNumber})
       expConfigData2{i,resultsXLSFNColNumber} = 'empty'; 
    end
end
[idx, label] = grp2idx(expConfigData2(:,resultsXLSFNColNumber));
[counts,labelN] = hist(idx,unique(idx));
% number of tasks * number of datasets = 28
valueCompleteExp = 28;
k = find(counts == valueCompleteExp);
% Unique group name (ResultsXLSFN) of completed experiments
groupExp = label(k);

%%
% Load ResultsControl xls file
% Insert alert that this script uses vars from anothers scripts
% workFolder
resultsControlPath = [workFolder '/' 'ResultsControl3.xlsx'];
% Load first sheet (Overview)
[~,~,rawRCOverview] = xlsread(resultsControlPath,1);
rcOverviewColHeadings = rawRCOverview(1,:);
rcOverview = rawRCOverview(2:end,:);
% Load second sheet (Detail)
[~,~,rawRCDetail] = xlsread(resultsControlPath,2);
rcDetailColHeadings = rawRCDetail(1,:);
rcDetail = cell2mat(rawRCDetail(2:end,:));

%%
% Copy proper data from ExperimentsConfig and Results xls files and put it
% in ResultsControl xls file
resultsPath = [workFolder '/Results/'];
for i = 1:length(groupExp)
    detailTempSheet = zeros(valueCompleteExp, length(rcDetailColHeadings));
    overviewTempSheet = cell(1,length(rcOverviewColHeadings));
    % Pick row indexes of the experiment group
    kk = find(strcmpi(groupExp(i), expConfigData2(:,resultsXLSFNColNumber)));
    
    if ~isempty(rcOverview), id = rcOverview{end,1} + 1;
    else id = i;
    end
    
    % Saves result id
    overviewTempSheet{1} = id;
    % Saves DR alg
    overviewTempSheet{2} = expConfigData2{kk(1),4};
    % Saves XLS file name
    overviewTempSheet{3} = expConfigData2{kk(1),resultsXLSFNColNumber};
    % Load group experiment xls file (ResultsMar2017..) - Summary sheet
    % For Cross-Validation
    [~,~,resultsSumRange] = xlsread([resultsPath overviewTempSheet{1,3}],...
        1,'Q73:Y74');
    % Cross validation general mean of SH(SH) and std
    overviewTempSheet{4} = resultsSumRange{1,1};
    overviewTempSheet{5} = resultsSumRange{2,1};
    % Cross validation general mean of SPD(SPD) and std
    overviewTempSheet{6} = resultsSumRange{1,5};
    overviewTempSheet{7} = resultsSumRange{2,5};
    % Cross validation general mean of SDBS(SDBS) and std
    overviewTempSheet{8} = resultsSumRange{1,9};
    overviewTempSheet{9} = resultsSumRange{2,9};
    % Load group experiment xls file (ResultsMar2017..) - Summary sheet
    % For Test set
    [~,~,resultsSumRange] = xlsread([resultsPath overviewTempSheet{1,3}],...
        1,'C73:K74');
    % Test general mean of SH(SH) and std
    overviewTempSheet{10} = resultsSumRange{1,1};
    overviewTempSheet{11} = resultsSumRange{2,1};
    % Test general mean of SPD(SPD) and std
    overviewTempSheet{12} = resultsSumRange{1,5};
    overviewTempSheet{13} = resultsSumRange{2,5};
    % Test general mean of SDBS(SDBS) and std
    overviewTempSheet{14} = resultsSumRange{1,9};
    overviewTempSheet{15} = resultsSumRange{2,9};
    for ii = 1:length(kk)
        % Save Results control ID
        detailTempSheet(ii,1) = id;
        % Save ExperimentID
        detailTempSheet(ii,2) = expConfigData2{kk(ii),1};
    end
    % Update loaded ResultsControl xls file var and
    % also update xls file
    nextOIndex = length(rcOverview(:,1)) + 2;
    nextOIndex = ['A' int2str(nextOIndex)];
    nextDIndex = size(rcDetail, 1) + 2;
    nextDIndex = ['A' int2str(nextDIndex)];
    
    xlwrite(resultsControlPath, overviewTempSheet, 1, nextOIndex);
    xlwrite(resultsControlPath, detailTempSheet, 2, nextDIndex);
    
    rcOverview = [rcOverview;overviewTempSheet];
    rcDetail = [rcDetail;detailTempSheet];
end
