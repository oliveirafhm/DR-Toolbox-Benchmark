% Trial setup generator using xls file
% Should be runned before everything to create various scenarios to be
% evaluated
% By F?bio Henrique (oliveirafhm@gmail.com)
% DR algorithm evaluation using machine learning
% 04/2017
%% Load scrips to write in xls file
clearvars;
poiPath = 'Third party codes/xlswrite/20130227_xlwrite/poi_library/';
javaaddpath([poiPath 'poi-3.8-20120326.jar']);
javaaddpath([poiPath 'poi-ooxml-3.8-20120326.jar']);
javaaddpath([poiPath 'poi-ooxml-schemas-3.8-20120326.jar']);
javaaddpath([poiPath 'xmlbeans-2.3.0.jar']);
javaaddpath([poiPath 'dom4j-1.6.1.jar']);
javaaddpath([poiPath 'stax-api-1.0.1.jar']);

%% Parameters initialization

% Load experiment config file
[configFileName, configPathName] = uigetfile('.xlsx', ...
    'Select (empty) experiment config file');
configFilePath = strcat(configPathName, configFileName);
[~,~,rawConfig] = xlsread(configFilePath,1);
% Change this variable every time that experiment configuration file has 
% changed (when added or removed a result column)
nResultsCols = 10;
expColHeadings = rawConfig(1,:);
expConfigData = cell2struct(rawConfig(2:end,1:end-nResultsCols),...
    expColHeadings(1:end-nResultsCols), 2);

nOfParams = length(expColHeadings(1:end-nResultsCols));

% Experiment ID
if ~isempty(expConfigData), id = expConfigData(end).ID + 1;
else id = 1;
end

[~,~,paramRange{1}] = xlsread(configFilePath,3); % PCA
[~,~,paramRange{2}] = xlsread(configFilePath,4); % Sammon
[~,~,paramRange{3}] = xlsread(configFilePath,5); % t-SNE

% paramFields = paramRange{1}(2:end,1);
% for i=1:3
% newExpConfigData(i) = cell2struct(paramRange{i}(2:end,2),...
%     paramFields);
% end

% First step of experiment config parse
nArrangements = [];
arrangement = 1;
for i = 1:length(paramRange)
    for j = 3:length(paramRange{i})
        if sum(paramRange{i}{j,2} == -1)==1 ...
                || sum(paramRange{i}{j,3} == 0)==1 ...
                || ischar(paramRange{i}{j,2})
            newExpConfigData1{i,j-1} = paramRange{i}{j,2};
        else
            newExpConfigData1{i,j-1} = ...
                paramRange{i}{j,2}:paramRange{i}{j,3}:paramRange{i}{j,4};
            arrangement = arrangement * length(newExpConfigData1{i,j-1});
        end
    end
    nArrangements = [nArrangements;arrangement];
    arrangement = 1;
end
%% Second step of experiment config parse

% Index of each array field from each DR alg
for i = 1:length(newExpConfigData1(:,1))
    kk{i} = find(cell2mat(paramRange{i}(2:end,5))==1)';
end

% Pick array fields and organize into cell structure
arrayFieldsCell = {};
arrayFieldsCellTemp = {};
for i = 1:length(paramRange)
    for j = 1:length(kk{i})
        arrayFieldsCellTemp{end+1,1} = paramRange{i}{kk{i}(j)+1,1};
        arrayFieldsCellTemp{end,2} = newExpConfigData1{i,kk{i}(j)};
    end
    arrayFieldsCell{i} = arrayFieldsCellTemp;
    arrayFieldsCellTemp = {};
end
% Cell + struct (same data as arrayFieldsCell)
arrayFieldsStruct = {};
for i = 1:length(arrayFieldsCell)
    for j = 1:length(arrayFieldsCell{i})
        field = arrayFieldsCell{i}{j,1};
        value = arrayFieldsCell{i}{j,2};
        arrayFieldsStruct{i}{j} = struct(field,value,'status',0,'col',...
            kk{i}(j));
    end
end

% Loop in arrayFieldsStruct to do all possible combinations (DR alg params)
for i = 1:length(arrayFieldsStruct)
    param = {};
    for j = 1:length(arrayFieldsStruct{i})
        fn = fieldnames(arrayFieldsStruct{i}{j});
        param{j} = arrayFieldsStruct{i}{j}.(fn{1});
    end
    % Alert: You should modify the source of combvec function before run it here
    allPossibilities{i} = combvec(param)';
end

% newExpConfigData = {sum(nArrangements),nOfParams};
newExpConfigData = {};
for i = 1:length(newExpConfigData1(:,1))
    newExpConfigData2 = {};
    for j = 1:nArrangements(i)        
        ii = 1:length(newExpConfigData1(1,:));
        % Get static values to compose partial experiment config table
        for k = ii(~ismember(ii,kk{i}))
            newExpConfigData2{j,k} = newExpConfigData1{i,k};
        end
        % Set Experiment ID        
        newExpConfigData2{j,1} = id;
        id = id + 1;
        % Get did arrangements an put in the partial experiment config
        % table
        for k = 1:length(kk{i})
          newExpConfigData2{j,kk{i}(k)} = allPossibilities{i}(j,k);
        end        
    end
    newExpConfigData = [newExpConfigData ; newExpConfigData2];
end

%% Save experiments config data in XLS file

% TODO: Fix whole section to allow append in a exist xls file
expFileName = ['ExperimentsConfig_' date '.xlsx'];
emptyExperimentsFilePath = [configPathName 'ExperimentsConfig_Empty.xlsx'];
copyfile(emptyExperimentsFilePath, [configPathName expFileName]);

xlwrite([configPathName expFileName], newExpConfigData, 1, 'A2');

%% End
clearvars;

