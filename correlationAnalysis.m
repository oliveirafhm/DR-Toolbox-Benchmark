% Author: Fabio Henrique (oliveirafhm@gmail.com)
% 04/2018
% DR parameters correlation analysis 

%% Load data
% Load parameter settings xls

%%
corrData_sammons = sammons;
% Pick only columns to be correlated
varsColumns = {'Iter';'LR';'DRErr';'CVAcc'};
corrDataSTable = array2table(corrData_sammons,...
    'VariableNames',varsColumns);
% Plot correlation data
baseFontSize = 16;
% Pearson's linear correlation coefficient
[RS, PValueS] = corrplot(corrDataSTable,'testR','on');
% title(['Sammons mapping parameters correlation'],'FontSize',baseFontSize);
% set(gca,'fontsize',baseFontSize+2);
% set(findobj(gca,'type','text'),'fontsize',baseFontSize);
% hold off; 

%%
corrData_tsne = t_sne;
% Pick only columns to be correlated
varsColumns = {'Iter';'LR';'Perp';'DRErr';'CVAcc'};
corrDataSTable = array2table(corrData_tsne,...
    'VariableNames',varsColumns);
% Plot correlation data
baseFontSize = 16;
% Pearson's linear correlation coefficient
[RT, PValueT] = corrplot(corrDataSTable,'testR','on');