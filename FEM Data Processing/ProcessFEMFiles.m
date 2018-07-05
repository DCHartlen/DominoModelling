%% Process FEM files
% 
% Script to process the GLSTAT file from domino simulations of dominoes.
% Includes final cropping, detrending, peak finding, and basic frequency analysis.
%
% Created by:  D.C. Hartlen, EIT
% Date:        ??-Mar-2018
% Modified by: D.C. Hartlen, EIT 
% Date:        18-Apr-2018

close all
clear
clc

screenSize = get( groot, 'Screensize' );

% Load Data
% [yy,Fs] = audioread('Domino_7pt5_1.m4a');
fileName = uigetfile({'*.*'},'Select GLSTAT File');
glstatData = glstat_parser(fileName);


xx = glstatData.Time;
yy = glstatData.ZVel;
Fs = (xx(2)-xx(1))^(-1);

figure('Name', 'Crop Full Data',...
    'OuterPosition',[0 0 screenSize(3) screenSize(4)])
plot(xx,yy)
yl = ylim;
ylim(yl)
xlabel('Time (s)')
ylabel('Amp')
title('Crop Data')
hold on

% Enable interactive cropping
[bx(1),~] = ginput(1);
line([bx(1),bx(1)],yl,'Color','k')
[bx(2),~] = ginput(1);
line([bx(2),bx(2)],yl,'Color','k')

bx = sort(bx);  % Sort bx to place smaller index first

% Crop data
xxCropped = xx(bx(1)*Fs:bx(2)*Fs);
% xx = xx-xx(1);
yyCropped = yy(bx(1)*Fs:bx(2)*Fs);

yyDetrend = PolyDetrend(xxCropped,yyCropped,3);

figure('Name', 'Processed Data',...
    'OuterPosition',[0 0 screenSize(3) screenSize(4)])
subplot(2,2,1)
hold on 
plot(xxCropped, yyCropped)
plot(xxCropped, yyDetrend)
% ylim([-1,1])
xlabel('Time (s)')
ylabel('Amp')
title('Cropped Data')
legend('Cropped','Cropped & Detrended')
 
% % Rectify data
% yy = yy.^2;
% 
%Develop and apply low pass filter to data
Fpass = 75;
Fstop = 300;
Apass = 1.0;
Astop = 65;
d = designfilt('lowpassiir', ...
  'PassbandFrequency',Fpass,'StopbandFrequency',Fstop, ...
  'PassbandRipple',Apass,'StopbandAttenuation',Astop, ...
  'DesignMethod','butter','SampleRate',Fs);
% fvtool(d)

yyDetrend = filtfilt(d,yyDetrend);
% 
% Plot Filtered Data
% yyDetrend = -yyDetrend;
subplot(2,2,3)
hold on
plot(xxCropped,smooth(yyDetrend))
xlabel('Time (s)')
ylabel('Amp')
title('Isolated Peaks')

% Find Peaks in filtered data
peakHeightThreshold = 5.0e-4;
peakSeperationThreshold = 0.030;
[peakVal,peakLoc,w,prom] = findpeaks(smooth(yyDetrend),Fs,...
                              'MinPeakDistance',peakSeperationThreshold,...
                              'MinPeakProminence', peakHeightThreshold);
hold on
plot(peakLoc+xxCropped(1),peakVal,'ro')
legend('Filtered Data','Peaks')

% Time derivative of impact time (peaks)
for i=1:length(peakLoc)-1
    delLoc(i) = peakLoc(i+1)-peakLoc(i);
end

% Plot Frequency
subplot(2,2,4)
plot(peakLoc(2:end),1./delLoc,'b-*')
hold on
meanFreq = mean(1./delLoc);
xl = xlim;
line(xl,meanFreq.*[1,1],'Color','r')
xlabel('Time (s)')
ylabel('Frequency (1/s)')
title(['Impact Frequency wrt Time: Mean Freq = ',num2str(meanFreq)])
% 
% % % Perform PSD on filtered data
% % [pxx,f] = periodogram(yyFilt, [],[],Fs);
% % figure()
% % plot(f,pxx)
% % xlim([-inf,100])
% 
% Perform PSD on unfiltered Data
[pxx,f] = periodogram(yyDetrend,hanning(length(yyDetrend)),[],Fs);
subplot(2,2,2)
hold on
plot(f(5:end),pxx(5:end))
xlabel('Frequency (1/s)')
ylabel('Amp')
title('PSD of cropped data')
xlim([0,inf])
% xlim([0,150])

