% Main file
% Uses PCAAnalysis and Classification m files
% By F?bio Henrique (oliveirafhm@gmail.com)
% 04/2017
%%
poiPath = 'Third party codes/xlswrite/20130227_xlwrite/poi_library/';
javaaddpath([poiPath 'poi-3.8-20120326.jar']);
javaaddpath([poiPath 'poi-ooxml-3.8-20120326.jar']);
javaaddpath([poiPath 'poi-ooxml-schemas-3.8-20120326.jar']);
javaaddpath([poiPath 'xmlbeans-2.3.0.jar']);
javaaddpath([poiPath 'dom4j-1.6.1.jar']);
javaaddpath([poiPath 'stax-api-1.0.1.jar']);
addpath(genpath('Third party codes'));

clc;
useConfigData = 1;
% Load experiment config file
if useConfigData == 1
    [configFileName, configPathName] = uigetfile('.xlsx', ...
        'Select experiment config file');
    configFilePath = strcat(configPathName, configFileName);
    [~,~,rawConfig] = xlsread(configFilePath,1);
    [~,~,colHeadings] = xlsread(configFilePath,2);
    colHeadingsTxt = colHeadings(1,:);
    colHeadingsLetters = colHeadings(2,:);
    % Change this variable every time that experiment configuration file has 
    % changed (when added or removed a result column)
    nResultsCols = 10; 
    expColHeadings = rawConfig(1,:);
    expConfigData = cell2struct(rawConfig(2:end,1:end-nResultsCols),...
        expColHeadings(1:end-nResultsCols), 2);
    % Loop through all config data
    for mainCounter = 1:length(expConfigData)
        if expConfigData(mainCounter).Status == 0
            tic
            % Save experiment datetime in ExperimentConfig xls file
            datetimeCol = colHeadingsLetters{ismember(colHeadingsTxt,'DateTime')};
            xlwrite(configFilePath, {datestr(datetime)}, 1, [datetimeCol int2str(mainCounter+1)]);
            PCAAnalysis;
            Classification;
            clearvars -except expConfigData mainCounter configFilePath ...
                useConfigData workFolder saveData colHeadingsLetters ...
                colHeadingsTxt;            
            expConfigData(mainCounter).Status = 1;
            % Update experiment status (to done) in ExperimentsConfig xls file
            statusCol = colHeadingsLetters{ismember(colHeadingsTxt,'Status')};
            xlwrite(configFilePath, 1, 1, [statusCol int2str(mainCounter+1)]);
            toc
        end
    end
    % Fit experiments data in ResultsControl xls file
    if saveData == 'y'
        ResultsControl;
    end
elseif useConfigData == 0
    PCAAnalysis;
    Classification;
end