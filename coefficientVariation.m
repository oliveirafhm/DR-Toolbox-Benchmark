% Compute Coefficient of Variation (CV)
% https://en.wikipedia.org/wiki/Coefficient_of_variation
% https://www.mathworks.com/matlabcentral/fileexchange/62972-getcv-x?focused=7658276&tab=function
% In response of:
% ‘The rest of samples, which is 10% as described in section 2.6, compose the test set.’
% How similar the same person was to himself/herself?
% Author: Fabio Henrique (oliveirafhm@gmail.com) - 16/09/2018

addpath(genpath('Third party codes'));

%% Load dataset
% Check code from PCAnalysis.m
load 'FeatMatrix-9-10-11.mat';

%% Creates variables to filter data and calculate CV
nGroups = 3;
nTasks = 4;
% Max number of samples of a specific class and task
maxN = 80 / nTasks;
z = cell(nGroups,1);
% Calculates CV
% Subject mean CVvalue x Task x Group
cvMean = NaN(maxN, nTasks, nGroups);
% Subject std CVvalue x Task x Group
cvStd = NaN(maxN, nTasks, nGroups);
% Do the job
for task = 1:4 % 1 - 4
    for i=1:nGroups
        if i == 1
            % Pick only non parkinson subjects
            iia{i} = find(FEATMATRIX(:,4) == 0 & FEATMATRIX(:,2) == task);
            t{i} = 'S_{H}';
            z{i} = FEATMATRIX(iia{i},6:end);
        elseif i == 2
            % Pick only parkinson subjects
            iia{i} = find(FEATMATRIX(:,4) == 1 & FEATMATRIX(:,5) == 0 ...
                & FEATMATRIX(:,2) == task);
            t{i} = 'S_{PD}';
            z{i} = FEATMATRIX(iia{i},6:end);
        elseif i == 3
            % Pick only dbs subjects
            iia{i} = find(FEATMATRIX(:,4) == 1 & FEATMATRIX(:,5) == 1 ...
                & FEATMATRIX(:,2) == task);
            t{i} = 'S_{DBS}';
            z{i} = FEATMATRIX(iia{i},6:end);
        end
        
        nSubjects = unique(FEATMATRIX(iia{i}, 1));
        subjectCounter = 1;
        for s = 1:length(nSubjects)
            iis = find(FEATMATRIX(iia{i},1) == nSubjects(s));
            cvValues = [];
            for nFeat = 1:length(z{i})
                cvValues(nFeat) = getCV(z{i}(iis, nFeat));
            end
            meanCV = mean(cvValues);
            stdCV = std(cvValues);
            cvMean(s, task, i) = meanCV;
            cvStd(s, task, i) = stdCV;
        end
    end
end
% Per task
grandMeanTaskGroup = nanmean(cvMean);
grandStdTaskGroup = nanmean(cvStd);
% openvar(grandMeanTaskGroup);
% openvar(grandStdTaskGroup);
% CV per group only
ggMeanCVGroup = mean(grandMeanTaskGroup);
ggStdCVGroup = mean(grandStdTaskGroup);
openvar('ggMeanCVGroup');
openvar('ggStdCVGroup');
