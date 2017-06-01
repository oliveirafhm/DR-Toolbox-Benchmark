% By F?bio Henrique (oliveirafhm@gmail.com)
% Run at the end of everything
% ...
% 06/2017

%% Line chart to compare overall results of each method (PCA, Sammon, t-SNE)

[resultsControlFileName, resultsControlPathName] = uigetfile('.xlsx', ...
    'Select results control xlsx file');
resultsControlFilePath = strcat(resultsControlPathName, resultsControlFileName);

[~,~,resultsControl] = xlsread(resultsControlFilePath,1);
resultsControlHeadings = resultsControl(1,1:end);
resultsControlData = resultsControl(2:end,1:end);