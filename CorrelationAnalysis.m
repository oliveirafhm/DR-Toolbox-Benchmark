%% Load data
% Load parameter settings xls

%%
corrDataS = sammons;
% Pick only columns to be correlated
varsColumns = {'Iter';'LR';'DRErr';'CVAcc'};
corrDataSTable = array2table(corrDataS,...
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
corrDataS = t_sne;
% Pick only columns to be correlated
varsColumns = {'Iter';'LR';'Perp';'DRErr';'CVAcc'};
corrDataSTable = array2table(corrDataS,...
    'VariableNames',varsColumns);
% Plot correlation data
baseFontSize = 16;
% Pearson's linear correlation coefficient
[RT, PValueT] = corrplot(corrDataSTable,'testR','on');