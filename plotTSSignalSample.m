% Author: Fabio Henrique (oliveirafhm@gmail.com)
% 08/04/2018
% Script to load a tremsen file from Alessandro's PhD study

%% Parte 1

addpath(genpath('adrianos_code'));
% Loading the LTFAT toolbox
% ltfatstart;

fs = 50;

% Load Annotation xls by hand
%[DATA]  = xlsread('Annotation.xlsx',  'DATA', 'A2:H153');
DATA = AnnotationS2;

subjects = unique(DATA(:,5));

%filename = ['Patient' num2str(subjects(s)) '.txt'];
%A = importdata([filename]); %carregando arquivo de dados
% Load tremsen file
[tsFileName, tsPathName] = uigetfile('.txt', ...
    'Select tremsen signal file');
filename = strcat(tsPathName, tsFileName);
A = importdata([filename]);

time = A.data(:,1);
s = str2num(tsFileName(8:9)); % subject id, pick it from file name

ii = find(DATA(:,5) == subjects(s));
hardwareType = unique(DATA(ii,7));

%% Sinal filtrado, normalizado e devidamento centralizado em
%torno do zero
[R] = EstimateResultant(hardwareType, A);
Nsinais = size(R,2);
if hardwareType == 1, pulseCol = 22; else pulseCol = 40; end 
kk = [1, 17, 25, 33]; % Sensor col
for k = kk
    %     [R(:,k), ~, ~] = GetInstaneousAtributes(R(:,k),fs);
    [a,~,~,~, ~, ~] = decomposeSig(R(:,k)); %decomposi??o do sinal original
    R(:,k) = R(:,k) - a;
    
    figure;
    plot(time, R(:,k), time, A.data(:,pulseCol));
end

%% Parte 2

% path = 'C:\Users\Adriano\Dropbox\Clinic\Patient 2/';
% filename = 'OFF A4 G2000 M2.txt';
[filename, path] = uigetfile('.txt', ...
    'Select file');
data = importdata([path  filename]); %load file
% data = data.data;

time = data.data(:,1); %time
marker1 = data.data(:,22);
marker2 = data.data(:,23);

%%
hardwareType = 1;
[R] = EstimateResultant(hardwareType, data);
fs = 50;
kk = [1, 17, 25, 33]; % Sensor col
for k = kk
    [R(:,k), ~, ~] = GetInstaneousAtributes(R(:,k),fs);    
%     figure;
%     plot(time, R(:,k), time, A.data(:,pulseCol));
end

figure;
%data from gyroscope
subplot(4,1,1);
plot(time,R(:,kk(1)));

%data from accelerometer
subplot(4,1,2);
plot(time,R(:,kk(2)));

%data from magnetometer
subplot(4,1,3);
plot(time,R(:,kk(3)));

%data from EMG 1
subplot(4,1,4);
plot(time,R(:,kk(4)));

%% Not used
figure;
%data from gyroscope
subplot(3,3,1);
plot(time,data(:,2)); 
hold on;
plot(time,data(:,3),'r'); 
plot(time,data(:,4),'g'); legend('X','Y','Z');
plot(time,marker1,'k');
plot(time,max(data(:,2))*(marker1),'y');
title('gyroscope 1');
xlim([min(time) max(time)]);

%data from accelerometer
subplot(3,3,4);
plot(time,data(:,5)); 
hold on;
plot(time,data(:,6),'r'); 
plot(time,data(:,7),'g'); legend('X','Y','Z');
plot(time,marker1,'k');
plot(time,max(data(:,7))*(marker1),'y');
title('accelerometer 1');
xlim([min(time) max(time)]);

%data from magnetometer
subplot(3,3,7);
plot(time,data(:,8)); 
hold on;
plot(time,data(:,9),'r'); 
plot(time,data(:,10),'g'); legend('X','Y','Z');
plot(time,marker1,'k');
plot(time,max(data(:,10))*(marker1),'y');
title('magnetometer 1');
xlim([min(time) max(time)]);

%data from gyroscope
subplot(3,3,2);
plot(time,data(:,11)); 
hold on;
plot(time,data(:,12),'r'); 
plot(time,data(:,13),'g'); legend('X','Y','Z');
plot(time,marker1,'k');
plot(time,max(data(:,13))*(marker1),'y');
title('gyroscope 2');
xlim([min(time) max(time)]);

%data from accelerometer
subplot(3,3,5);
plot(time,data(:,14)); 
hold on;
plot(time,data(:,15),'r'); 
plot(time,data(:,16),'g'); legend('X','Y','Z');
plot(time,marker1,'k');
plot(time,max(data(:,16))*(marker1),'y');
title('accelerometer 2');
xlim([min(time) max(time)]);

%data from magnetometer
subplot(3,3,8);
plot(time,data(:,17)); 
hold on;
plot(time,data(:,18),'r'); 
plot(time,data(:,19),'g'); legend('X','Y','Z');
plot(time,marker1,'k');
plot(time,max(data(:,19))*(marker1),'y');
title('magnetometer 2');
xlim([min(time) max(time)]);

%data from EMG 1
subplot(3,3,3);
plot(time,data(:,20)); 
hold on;
plot(time,data(:,21),'r'); 
legend('EMG1','EMG2');
plot(time,marker1,'k');
plot(time,max(data(:,21))*(marker1),'y');
title('EMG');
xlim([min(time) max(time)]);
