%% Plot decision boundary in scatter plot
% Load .mat data before proceed (_OosModelData.mat)
addpath(genpath('Third party codes'));
[matFileName, matPathName] = uigetfile('/Users/oliveirafhm/Dropbox/Compartilhada Alessando e Fabio/Paper Neural Computing and Applications 2016-2/Paper doc word/images/data_visualization/*.mat', ...
    'Select data file related to scatter plot');
matFilePath = strcat(matPathName, matFileName);
load(matFilePath);

pcaIndex = strfind(matFileName,'PCA');
tsneIndex = strfind(matFileName,'t-SNE');
if pcaIndex > 0
    [trainDataPCA, mappingInfoPCA] = compute_mapping(trainData, 'PCA', 2);
    P = trainDataPCA';
else
    [trainDataPCA, mappingInfoPCA] = compute_mapping(trainData, 'PCA', 0.9);
    P = bestTrainedN.net(trainDataPCA');
end

T = trainLabels';

%% create a neural network and train

ffnet = feedforwardnet([4 3]);
ffnet.trainParam.epochs = 2000;
ffnet.divideParam.trainRatio = 1;
ffnet.divideParam.valRatio   = 0;
ffnet.divideParam.testRatio  = 0;
[ffnet,tr,Y,E] = train(ffnet,P,T);
% view(ffnet)

%% evaluate performance: decoding network response
% [m,i] = max(T); % target class
% [m,j] = max(Y); % predicted class
% N = length(Y);  % number of all samples
% k = 0;          % number of missclassified samples
% if find(i-j),   % if there exist missclassified samples
%     k = length(find(i-j)); % get a number of missclassified samples
% end
% fprintf('Correct classified samples: %.1f%% samples\n', 100*(N-k)/N);

%% generate a grid
% Load figure before proceed (_ScatterPlot.fig)
splitPoints = find(matFileName == '_');
figName = [matFileName(1:splitPoints(end)) 'ScatterPlot.fig'];
uiopen([matPathName figName],1);

f1 = figure(1);
xl = xlim; yl = ylim;
if tsneIndex > 0, spanStep = .05; else spanStep = .1; end
spanX = xl(1):spanStep:xl(2);
spanY = yl(1):spanStep:yl(2);
[P1,P2] = meshgrid(spanX,spanY);
pp = [P1(:) P2(:)]';

% simulate neural network on a grid
aa = ffnet(pp);
aa = round(aa);
%% plot classification regions based on MAX activation
hold on;
red = [1 0 0];
black = [0 0 0];
blue = [0 0 1];
map = [blue;black;red];
% m = mesh(P1,P2,reshape(double(ismember(aa,1)),length(spanY),length(spanX)));
%f2 = figure;
z = reshape(aa,length(spanY),length(spanX));
m = mesh(P1,P2,z,'EdgeAlpha',0.1,'FaceAlpha',0.1);
view(2);
colormap(f1,map);
%Returns handles to the patch and line objects
% chi=get(gca, 'Children');
%Reverse the stacking order so that the patch overlays the line
% set(gca, 'Children',flipud(chi));

%%
%decisionB = findobj(2, 'type', 'surface');
%copyobj(decisionB,findobj(1,'type','axes'));

%% Save figure

pbaspect([1.25 1 1]);
qstring = 'Would you like to save this figure?';
choice = questdlg(qstring,'Decision boundary map',...
    'Yes','No','Yes');

if strcmp(choice,'Yes')
    figurePathName = [matFileName(1:splitPoints(end)) ...
        'ScatterPlot_WithBoundary.fig'];
    savefig([matPathName figurePathName]);
end
close all;
% clearvars;
clear all;

