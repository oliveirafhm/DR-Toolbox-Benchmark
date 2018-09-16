% 1 - Put Clinic folder in matlab path
% 2 - Script used to plot sample signals from tremsen

%MacOSX
%rootFolder = '/Users/oliveirafhm/OneDrive/Doutorado/Compartilhada_Alessando_Fabio/';
%Windows
rootFolder = 'C:\Users\olive\OneDrive\Doutorado\Compartilhada_Alessando_Fabio\';
addpath(genpath(rootFolder));

%% Delete all variables and close all opened figures

clear all;
close all;

addpath(genpath('ltfat-2.1.0-win64'));

%% Loading the LTFAT toolbox
ltfatstart;

%% Load data file and select subject
% MacOSX
% folder = '/Users/oliveirafhm/OneDrive/Doutorado/Compartilhada_Alessando_Fabio/Clinic_09_03_2016/';
% Windows
folder = 'C:\Users\olive\OneDrive\Doutorado\Compartilhada_Alessando_Fabio\Clinic_09_03_2016\';

Subject = 18;
tSubject = ['Patient' num2str(Subject)];
% MacOSX
%A = importdata([folder tSubject '/' tSubject '.txt']);
%Windows
A = importdata([folder tSubject '\' tSubject '.txt']);

fs = 50; %sampling frequency in Hz

t = A.data(:,1); % time in seconds
t_bkp = t;

%% Load annotation
[Annotation]  = xlsread('Annotation.xlsx',  'DATA', 'A2:H761');
% column 1 = wo -> tempo inicial da janela (segundos)
% column 2 = wf -> tempo final da janela (segundos)
% column 3 = Task -> tarefa (1 = pin?a, 2 = m?o at? o nariz, 3 = rota??o da m?o, 4 =  m?o est?tica)
% column 4 = Stimulus intensity -> intensidade do estimulo ( 1 = 0%; 2 = 25%; 3 = 50%; 4 = 75%; 5 = 100%)
% column 5 = Subject -> identificador do sujeito (4 a 42)
% column 6 = Pathological -> 1 = presen?a de doen?a (tremor); 0 = aus?ncia de tremor
% column 7 = File format -> 1 = primeira vers?o do hardware; 2 = segunda vers?o do hardware; 3 = terceira vers?o do hardware
% column 8 = DBS -> 1 = presen?a de estimulador; 0 = aus?ncia

%%%%%%

%load 'Annotation.mat';
ii = find(Annotation(:,5) == Subject);
M = Annotation(ii,:);
c = {'c' 'k' 'g' 'b' 'r'}; %each colour identifies a particular stimulus intensity 


%% estimate resulting componenent for each sensor

if(M(1,7) == 1)
    
    gyroS1X = A.data(:,2);
    gyroS1Y = A.data(:,3);
    gyroS1Z = A.data(:,4);
    gyroS1 = sqrt(A.data(:,2).^2 + A.data(:,3).^2 + A.data(:,4).^2); 
    
    
    gyroS2X = A.data(:,5);
    gyroS2Y = A.data(:,6);
    gyroS2Z = A.data(:,7);
    gyroS2 = sqrt(A.data(:,5).^2 + A.data(:,6).^2 + A.data(:,7).^2); 
    
    accS1 = sqrt(A.data(:,8).^2 + A.data(:,9).^2 + A.data(:,10).^2); 
    accS2 = sqrt(A.data(:,11).^2 + A.data(:,12).^2 + A.data(:,13).^2); 

    magS1 = sqrt(A.data(:,14).^2 + A.data(:,15).^2 + A.data(:,16).^2); 
    magS2 = sqrt(A.data(:,17).^2 + A.data(:,18).^2 + A.data(:,19).^2); 

    EMG1 = A.data(:,20);
    EMG2 = A.data(:,21);
    
    pulsoA = A.data(:,22); % pulseA
    
end
    
if(M(1,7) == 2)
        
    gyroS1X = A.data(:,2);
    gyroS1Y = A.data(:,3);
    gyroS1Z = A.data(:,4);
    gyroS1 = sqrt(A.data(:,2).^2 + A.data(:,3).^2 + A.data(:,4).^2); 
    
    
    gyroS2X = A.data(:,8);
    gyroS2Y = A.data(:,9);
    gyroS2Z = A.data(:,10);
    gyroS2 = sqrt(A.data(:,8).^2 + A.data(:,9).^2 + A.data(:,10).^2); 
    
    accS1 = sqrt(A.data(:,14).^2 + A.data(:,15).^2 + A.data(:,16).^2); 
    accS2 = sqrt(A.data(:,20).^2 + A.data(:,21).^2 + A.data(:,22).^2); 

    magS1 = sqrt(A.data(:,26).^2 + A.data(:,27).^2 + A.data(:,28).^2); 
    magS2 = sqrt(A.data(:,32).^2 + A.data(:,33).^2 + A.data(:,34).^2); 

    EMG1 = A.data(:,38);
    EMG2 = A.data(:,39);
    
    pulsoA = A.data(:,40); % pulseA

end


if(M(1,7) == 3)
        
    gyroS1X = A.data(:,2);
    gyroS1Y = A.data(:,3);
    gyroS1Z = A.data(:,4);
    gyroS1 = sqrt(A.data(:,2).^2 + A.data(:,3).^2 + A.data(:,4).^2); 
    
    
    gyroS2X = A.data(:,5);
    gyroS2Y = A.data(:,6);
    gyroS2Z = A.data(:,7);
    gyroS2 = sqrt(A.data(:,5).^2 + A.data(:,6).^2 + A.data(:,7).^2); 
    
    accS1 = sqrt(A.data(:,14).^2 + A.data(:,15).^2 + A.data(:,16).^2); 
    accS2 = sqrt(A.data(:,17).^2 + A.data(:,18).^2 + A.data(:,19).^2); 

    magS1 = sqrt(A.data(:,26).^2 + A.data(:,27).^2 + A.data(:,28).^2); 
    magS2 = sqrt(A.data(:,29).^2 + A.data(:,30).^2 + A.data(:,31).^2); 

    EMG1 = A.data(:,38);
    EMG2 = A.data(:,39);
    
    pulsoA = A.data(:,40); % pulseA

end


%% Plot raw data
figure; 
[h1] = subplot(4,1,1);
plot(t,gyroS1);
hold on; plot(t,gyroS2,'r');
plot(t,max(gyroS2)*pulsoA/max(pulsoA),'k');
title('Gyroscope'); 
ylabel('degrees/s');

for i=1:1:size(M,1),
     plot([M(i,1) M(i,2)], [max(gyroS2) max(gyroS2)],c{M(i,4)}); 
end

xlim([t(1) t(end)]);

[h2] = subplot(4,1,2);
plot(t,accS1);
hold on; plot(t,accS2,'r');
plot(t,max(accS2)*pulsoA/max(pulsoA),'k');
title('Accelerometer');
ylabel('g');
for i=1:1:size(M,1),
     plot([M(i,1) M(i,2)], [max(accS2) max(accS2)],c{M(i,4)});
end
xlim([t(1) t(end)]);

[h3] = subplot(4,1,3);
plot(t,magS1);
hold on; plot(t,magS2,'r');
plot(t,max(magS2)*pulsoA/max(pulsoA),'k');
for i=1:1:size(M,1),
     plot([M(i,1) M(i,2)], [max(magS2) max(magS2)],c{M(i,4)});
end
xlim([t(1) t(end)]);
title('Magnetometer');
ylabel('gauss');

[h4] = subplot(4,1,4);
plot(t,EMG1);
hold on; plot(t,EMG2,'r');
plot(t,max(EMG1)*pulsoA/max(pulsoA),'k');
for i=1:1:size(M,1),
     plot([M(i,1) M(i,2)], [max(EMG1) max(EMG1)],c{M(i,4)});
end
xlim([t(1) t(end)]);
title('EMG');
ylabel('V');
xlabel('time (s)');
linkaxes([h1,h2,h3,h4],'x') ;

%% Plot normalized and filtered data
figure; 

Ymin = 0; %escala m?nima para normaliza??o
Ymax = 1; %escala m?xima para normaliza??o

[h1] = subplot(4,1,1);
[gyroS1_] = FiltSig(gyroS1,fs); %filtered signal
[gyroS2_] = FiltSig(gyroS2,fs); %filtered signal


NormGyro1 = convScale(min(gyroS1_),max(gyroS1_),gyroS1_,Ymin,Ymax);
NormGyro2 = convScale(min(gyroS2_),max(gyroS2_),gyroS2_,Ymin,Ymax);
plot(t,NormGyro1,'b');
hold on; plot(t,NormGyro2,'r');

%[pulsoA_]= convScale(min(pulsoA),max(pulsoA),pulsoA,min(NormGyro2),max(NormGyro2));
[pulsoA_]= convScale(min(pulsoA),max(pulsoA),pulsoA,Ymin,Ymax);
plot(t,pulsoA_,'k');
for i=1:1:size(M,1),
     plot([M(i,1) M(i,2)], [max(pulsoA_) max(pulsoA_)],c{M(i,4)});
end
xlim([t(1) t(end)]); title('Gyroscope - Filtered and Normalized'); 
ylim([Ymin 1.3*Ymax]); 


[h2] = subplot(4,1,2);
[accS1_] = FiltSig(accS1,fs); %filtered signal
[accS2_] = FiltSig(accS2,fs); %filtered signal

NormAccS1 = convScale(min(accS1_),max(accS1_),accS1_,Ymin,Ymax);
NormAccS2 = convScale(min(accS2_),max(accS2_),accS2_,Ymin,Ymax);
plot(t,NormAccS1,'b');
hold on; plot(t,NormAccS2,'r');
[pulsoA_]= convScale(min(pulsoA),max(pulsoA),pulsoA,min(NormAccS2),max(NormAccS2));
plot(t,pulsoA_,'k');
for i=1:1:size(M,1),
     plot([M(i,1) M(i,2)], [max(pulsoA_) max(pulsoA_)],c{M(i,4)});
end
xlim([t(1) t(end)]);title('Accelerometer - Filtered and Normalized');
ylim([Ymin 1.3*Ymax]); 



[h3] = subplot(4,1,3);
[magS1_] = FiltSig(magS1,fs); %filtered signal
[magS2_] = FiltSig(magS2,fs); %filtered signal
NormMagS1 = convScale(min(magS1_),max(magS1_),magS1_,Ymin,Ymax);
NormMagS2 = convScale(min(magS2_),max(magS2_),magS2_,Ymin,Ymax);
plot(t,NormMagS1,'b');
hold on; plot(t,NormMagS2,'r');
[pulsoA_]= convScale(min(pulsoA),max(pulsoA),pulsoA,min(NormMagS2),max(NormMagS2));
plot(t,pulsoA_,'k');
for i=1:1:size(M,1),
     plot([M(i,1) M(i,2)], [max(pulsoA_) max(pulsoA_)],c{M(i,4)}');
end
xlim([t(1) t(end)]);title('Magnetometer - Filtered and Normalized');
ylim([Ymin 1.3*Ymax]); 


[h4] = subplot(4,1,4);
[EMG1_] = FiltSig(EMG1,fs); %filtered signal
[EMG2_] = FiltSig(EMG2,fs); %filtered signal

NormEMG1 = convScale(min(EMG1_),max(EMG1_),EMG1_,Ymin,Ymax);
NormEMG2 = convScale(min(EMG2_),max(EMG2_),EMG2_,Ymin,Ymax);
plot(t,NormEMG1,'b');
hold on; plot(t,NormEMG2,'r');
[pulsoA_]= convScale(min(pulsoA),max(pulsoA),pulsoA,min(NormEMG2),max(NormEMG2));
plot(t,pulsoA_,'k');
for i=1:1:size(M,1),
     plot([M(i,1) M(i,2)], [max(pulsoA_) max(pulsoA_)],c{M(i,4)});
end
xlim([t(1) t(end)]);title('EMG - Filtered and Normalized'); xlabel('time (s)');
ylim([Ymin 1.3*Ymax]); 

set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold');

set(findall(gcf,'type','axes'),'FontSize',20);


linkaxes([h1,h2,h3,h4],'xy') ;

%% Plot normalized and filtered data (first repetition)
wn = [0 46];
t = t_bkp(t_bkp>=wn(1) & t_bkp<=wn(2));
figure;

Ymin = 0; %escala m?nima para normaliza??o
Ymax = 1; %escala m?xima para normaliza??o

[h1] = subplot(4,1,1);
[gyroS1_] = FiltSig(gyroS1,fs); %filtered signal
[gyroS2_] = FiltSig(gyroS2,fs); %filtered signal

NormGyro1 = convScale(min(gyroS1_),max(gyroS1_),gyroS1_,Ymin,Ymax);
NormGyro2 = convScale(min(gyroS2_),max(gyroS2_),gyroS2_,Ymin,Ymax);
NormGyro1 = NormGyro1(t>=wn(1) & t<=wn(2));
NormGyro2 = NormGyro2(t>=wn(1) & t<=wn(2));
plot(t,NormGyro1,'b');
hold on; %plot(t,NormGyro2,'r');

%[pulsoA_]= convScale(min(pulsoA),max(pulsoA),pulsoA,min(NormGyro2),max(NormGyro2));
[pulsoA_]= convScale(min(pulsoA),max(pulsoA),pulsoA,Ymin,Ymax);
pulsoA_ = pulsoA_(t>=wn(1) & t<=wn(2));
plot(t,pulsoA_,'k');
% for i=1:1:size(M,1),
%      plot([M(i,1) M(i,2)], [max(pulsoA_) max(pulsoA_)],c{M(i,4)});
% end
xlim([t(1) t(end)]); title('Gyroscope - Filtered and Normalized'); 
ylim([Ymin 1.1*Ymax]);

[h2] = subplot(4,1,2);
[accS1_] = FiltSig(accS1,fs); %filtered signal
[accS2_] = FiltSig(accS2,fs); %filtered signal

NormAccS1 = convScale(min(accS1_),max(accS1_),accS1_,Ymin,Ymax);
NormAccS2 = convScale(min(accS2_),max(accS2_),accS2_,Ymin,Ymax);
NormAccS1 = NormAccS1(t>=wn(1) & t<=wn(2));
NormAccS2 = NormAccS2(t>=wn(1) & t<=wn(2));
plot(t,NormAccS1,'b');
hold on; %plot(t,NormAccS2,'r');
[pulsoA_]= convScale(min(pulsoA),max(pulsoA),pulsoA,min(NormAccS2),max(NormAccS2));
pulsoA_ = pulsoA_(t>=wn(1) & t<=wn(2));
plot(t,pulsoA_,'k');
% for i=1:1:size(M,1),
%      plot([M(i,1) M(i,2)], [max(pulsoA_) max(pulsoA_)],c{M(i,4)});
% end
xlim([t(1) t(end)]);title('Accelerometer - Filtered and Normalized');
ylim([Ymin 1.1*Ymax]); 

[h3] = subplot(4,1,3);
[magS1_] = FiltSig(magS1,fs); %filtered signal
[magS2_] = FiltSig(magS2,fs); %filtered signal
NormMagS1 = convScale(min(magS1_),max(magS1_),magS1_,Ymin,Ymax);
NormMagS2 = convScale(min(magS2_),max(magS2_),magS2_,Ymin,Ymax);
NormMagS1 = NormMagS1(t>=wn(1) & t<=wn(2));
NormMagS2 = NormMagS2(t>=wn(1) & t<=wn(2));
plot(t,NormMagS1,'b');
hold on; %plot(t,NormMagS2,'r');
[pulsoA_]= convScale(min(pulsoA),max(pulsoA),pulsoA,min(NormMagS2),max(NormMagS2));
pulsoA_ = pulsoA_(t>=wn(1) & t<=wn(2));
plot(t,pulsoA_,'k');
% for i=1:1:size(M,1),
%      plot([M(i,1) M(i,2)], [max(pulsoA_) max(pulsoA_)],c{M(i,4)}');
% end
xlim([t(1) t(end)]);title('Magnetometer - Filtered and Normalized');
ylim([Ymin 1.1*Ymax]);

% ----
% [h4] = subplot(4,1,4);
% plot(t,EMG1);
% hold on; plot(t,EMG2,'r');
% plot(t,max(EMG1)*pulsoA/max(pulsoA),'k');
% for i=1:1:size(M,1),
%      plot([M(i,1) M(i,2)], [max(EMG1) max(EMG1)],c{M(i,4)});
% end
% xlim([t(1) t(end)]);
% title('EMG');
% ylabel('V');
% xlabel('time (s)');
% ----

[h4] = subplot(4,1,4);
% [EMG1_] = FiltSig(EMG1,fs); %filtered signal
% [EMG2_] = FiltSig(EMG2,fs); %filtered signal
% 
% NormEMG1 = convScale(min(EMG1_),max(EMG1_),EMG1_,Ymin,Ymax);
% NormEMG2 = convScale(min(EMG2_),max(EMG2_),EMG2_,Ymin,Ymax);
NormEMG1 = EMG1;
NormEMG2 = EMG2;
NormEMG1 = NormEMG1(t>=wn(1) & t<=wn(2));
NormEMG2 = NormEMG2(t>=wn(1) & t<=wn(2));
plot(t,NormEMG1,'b');
hold on; %plot(t,NormEMG2,'r');
[pulsoA_]= convScale(min(pulsoA),max(pulsoA),pulsoA,min(NormEMG2),max(NormEMG2));
pulsoA_ = pulsoA_(t>=wn(1) & t<=wn(2));
plot(t,pulsoA_,'k');
% for i=1:1:size(M,1),
%      plot([M(i,1) M(i,2)], [max(pulsoA_) max(pulsoA_)],c{M(i,4)});
% end
xlim([t(1) t(end)]);title('EMG Envelope'); xlabel('time (s)');
ylabel('V');
% ylim([Ymin 1.1*Ymax]);

set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold');

set(findall(gcf,'type','axes'),'FontSize',20);

linkaxes([h1,h2,h3],'xy'); % ,h4

%% spectrum analysis
figure;
OrderAr = 20; %order of the AR model
%signal filtering
Lo = 1; %lower cutoff frequency in Hz
Lu = fs/2; %upper cutoff frequency in Hz
%order = 5; %order

% y = NormGyro1, NormGyro2, NormAccS1, NormAccS2, NormMagS1,NormMagS2, NormEMG1, NormEMG2

y = NormGyro1;

[h1] = subplot(3,1,1);
plot(t,y); hold on;
ylim([Ymin 1.3*Ymax]); 

plot(t,pulsoA_,'k');
for i=1:1:size(M,1),
     plot([M(i,1) M(i,2)], [max(pulsoA_) max(pulsoA_)],c{M(i,4)}');
end

xlim([t(1) t(end)]);
xlabel('time (s)');

if (M(1,6) == 1)
    
    pathol = 'Tremor + DBS';
    
else
    
    pathol = 'Healthy';
end
    
    
title(['Subject: ' int2str(Subject) ' (' pathol ')' ]);



[h2] = subplot(3,1,2);
% spectrogram(y,128,120,128,fs,'yaxis');
[s1,f1,t1,p1] = spectrogram((y-mean(y))/std(y),128,120,[0.1:0.1:fs/2],fs);


%yy = spline(t1,s1,xx);




image(t1,f1,abs(s1)); colormap hot; set(gca,'Ydir','Normal');


[q,nd] = max(10*log10(p1));

hold on
plot3(t1,f1(nd),q,'g','linewidth',4)
% hold off
% %colorbar
% view(2)

% colormap hot
% colorbar(h2,'off');
%sgram(y,fs,90,'wlen',round(500/1000*fs),'nocolorbar'); 
hold on;
xlabel('time (s)');
ylabel('frequency (Hz)');

subplot(3,1,3);
[Pxx,F] = pburg(y,OrderAr,[Lo:(Lu-Lo)/500:Lu],fs) ;
plot(F,Pxx/sum(Pxx));

ylabel('energy (normalized)');
xlabel('frequency (Hz)');

%linkaxes([h1,h2],'x') ;


[LocalMaxV,LocalMaxI, LocalMinV,LocalMinI]= extrema(Pxx/sum(Pxx));
hold on; plot(F(LocalMaxI),LocalMaxV,'r*');


set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold');

set(findall(gcf,'type','axes'),'FontSize',20);

%% Rastreamento da frequencia fundamental
x = magS1;

[diffFreq,theta_curve,Xn] = FreqTrack(x,fs);
%x = Xn;


figure;
h1 = subplot(3,1,1);

plot(t,Xn); hold on;

title('Normalized and filtered signal');
xlabel('time (s)');


for i=1:1:size(M,1),
     plot([M(i,1) M(i,2)], [max(Xn) max(Xn)],c{M(i,4)});
end
xlim([t(1) t(end)]);
ylim([min(Xn) 1.1*max(Xn)]);

h2 = subplot(3,1,2);

plot(t, theta_curve,'r');hold on;

for i=1:1:size(M,1),
     plot([M(i,1) M(i,2)], [max(theta_curve) max(theta_curve)],c{M(i,4)});
end
xlim([t(1) t(end)]);
ylim([min(theta_curve) 1.1*max(theta_curve)]);

title('Fundamental frequency - f_o(t)');
xlabel('time (s)');
ylabel('Hz');

linkaxes([h1 h2], 'x');



% Estimativa de oscila??es de alta frequencia

h3 = subplot(3,1,3);

plot(t(1:end-1), diffFreq,'r');hold on;

for i=1:1:size(M,1),
     plot([M(i,1) M(i,2)], [max(diffFreq) max(diffFreq)],c{M(i,4)});
end
xlim([t(1) t(end)]);
ylim([min(diffFreq) 1.1*max(diffFreq)]);


title('f_h(t) = f_o(t)/dt');
xlabel('time (s)');
ylabel('Hz/s');

linkaxes([h1 h2 h3], 'x');


set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold');

set(findall(gcf,'type','axes'),'FontSize',20);

figure;

[x2N]= convScale(min(x),max(x),x,0,1);

plot(t,x2N);

hold on;

[diffFreqN]= convScale(min(diffFreq),max(diffFreq),diffFreq,0,1);
plot(t(1:end-1),diffFreqN,'r')


%% Filtragem adequada do sinal original (somente o sinal original conforme especificado abaixo) e estimativa da frequencia instantanea

%  gyroS1 gyroS2 accS1 accS2 magS1 magS2 EMG1 EMG2
xdt = gyroS1;

%vel = cumtrapz(xdt) * 1/fs + 0;
[XdtFiltered, yi, amp] = GetInstaneousAtributes(xdt,fs);

figure;

%Sinal original
h1 = subplot(4,1,1);
plot(t,xdt); 
title('Raw signal');

hold on;
for i=1:1:size(M,1),
     plot([M(i,1) M(i,2)], [max(xdt) max(xdt)],c{M(i,4)});
end


%Sinal filtrado
h2 = subplot(4,1,2);
plot(h2,t,XdtFiltered,'b'); hold on;
plot(h2,[t(1) t(end)],[0 0], 'k');
title('Filtered and normalized signal');

for i=1:1:size(M,1),
     plot([M(i,1) M(i,2)], [max(XdtFiltered) max(XdtFiltered)],c{M(i,4)});
end


%Frequencia instantanea
h3 = subplot(4,1,3); 
plot(t,yi); hold on;
%plot(h3,[t(1) t(end)],[mean(yi) mean(yi)], 'k');


for i=1:1:size(M,1),
     plot([M(i,1) M(i,2)], [max(yi) max(yi)],c{M(i,4)});
end

ylim([0 fs/2]);

title('Instantaneous frequency');
ylabel('frequency (Hz)');


h4 = subplot(4,1,4);

plot(t,amp);hold on;
%plot(h4,[t(1) t(end)],[mean(amp) mean(amp)], 'k');


for i=1:1:size(M,1),
     plot([M(i,1) M(i,2)], [max(amp) max(amp)],c{M(i,4)});
end


title('Instantaneous amplitude');
xlabel('time (s)');

linkaxes([h1 h2 h3 h4],'x');

xlim([t(1) t(end)]);

set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold');

set(findall(gcf,'type','axes'),'FontSize',20);