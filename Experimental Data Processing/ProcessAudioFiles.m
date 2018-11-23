%% Process Audio Files
% 
% Script to process sound recordings of dominoes. Includes final cropping,
% filtering, peak finding, and basic frequency analysis.
%
% Created by:  D.C. Hartlen, EIT
% Date:        ??-Mar-2018
% Modified by: D.C. Hartlen, EIT 
% Date:        17-Aug-2018

close all
clear
clc

screenSize = get( groot, 'Screensize' );

% Load Data
% [yy,Fs] = audioread('Domino_7pt5_1.m4a');
[fileName, pathname, ~] = uigetfile({'*.*'},'Select Audio File');
[yy,Fs] = audioread([pathname fileName]);


xx = linspace(0,length(yy)/Fs,length(yy))';

figure('Name', 'Crop Full Data',...
    'OuterPosition',[0 0 screenSize(3) screenSize(4)])
plot(xx,yy)
ylim([-1,1])
xlabel('Time (s)')
ylabel('Amp')
title('Crop Data')
hold on

% Enable interactive cropping
[bx(1),~] = ginput(1);
line([bx(1),bx(1)],[-1.5,1.5],'Color','k')
[bx(2),~] = ginput(1);
line([bx(2),bx(2)],[-1.5,1.5],'Color','k')
bx = sort(bx);

% bx = [1.43864406779661,1.80745762711864]

% Crop data
xx = xx(bx(1)*Fs:bx(2)*Fs);
% xx = xx-xx(1);
yy = yy(bx(1)*Fs:bx(2)*Fs);

figure('Name', 'Analysed Data',...
    'OuterPosition',[0 0 screenSize(3) screenSize(4)])
subplot(2,2,1)
plot(xx,yy)
ylim([-1,1])
xlabel('Time (s)')
ylabel('Amp')
title('Crop Data')

% Rectify data
yy = yy.^2;

%Develop and apply low pass filter to data
Fpass = 100;
Fstop = 300;
Apass = 0.25;
Astop = 75;
d = designfilt('lowpassfir', ...
  'PassbandFrequency',Fpass,'StopbandFrequency',Fstop, ...
  'PassbandRipple',Apass,'StopbandAttenuation',Astop, ...
  'DesignMethod','equiripple','SampleRate',Fs);
% fvtool(d)

yyFilt = filtfilt(d,yy);

% Plot Filtered Data
subplot(2,2,3)
hold on
plot(xx,yyFilt)
xlabel('Time (s)')
ylabel('Amp')
title('Isolated Peaks')

% Find Peaks in filtered data
peakHeightThreshold = 0.0000
5;
peakSeperationThreshold = 0.006;
[peakVal,peakLoc,w,prom] = findpeaks(yyFilt,Fs,...
                              'MinPeakDistance',peakSeperationThreshold,...
                              'MinPeakProminence',peakHeightThreshold);
hold on
plot(peakLoc(1:min(32,end))+xx(1),peakVal(1:min(32,end)),'ro')
legend('Filtered Data','Peaks')


% Time derivative of impact time (peaks)
for i=1:length(peakLoc)-1
    delLoc(i) = peakLoc(i+1)-peakLoc(i);
end

% Plot Frequency
subplot(2,2,4)
plot(peakLoc(2:min(32,end)),1./delLoc(1:min(31,end)),'b-*')
hold on
meanFreq = mean(1./delLoc(1:min(31,end)));
xl = xlim;
line(xl,meanFreq.*[1,1],'Color','r')
xlabel('Time (s)')
ylabel('Frequency (1/s)')
title(['Impact Frequency wrt Time: Mean Freq = ',num2str(meanFreq)])

% % Perform PSD on filtered data
% [pxx,f] = periodogram(yyFilt, [],[],Fs);
% figure()
% plot(f,pxx)
% xlim([-inf,100])

% Perform PSD on unfiltered Data
[pxx,f] = periodogram(yy,rectwin(length(yy)),2*2^nextpow2(length(yy)),Fs);
subplot(2,2,2)
hold on
plot(f,pxx)
xlabel('Frequency (1/s)')
ylabel('Amp')
title('PSD of cropped data')
xlim([0,150])

