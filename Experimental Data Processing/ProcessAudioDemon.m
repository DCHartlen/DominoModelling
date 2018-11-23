%% Process Audio DEMON
% 
% Utilizes Detection of Envelope Modulation of Noise to capture domino
% impacts based on how the shape of the signal envelop changes, not from
% individual peaks in the audio signal.
%
% Created by:  D.C. Hartlen, EIT
% Date:        17-Aug-2018
% Modified by:  
% Date:        

close all
clear
clc

%% Crop data
screenSize = get( groot, 'Screensize' );

% Load Data
[fileName, pathname, ~] = uigetfile({'*.*'},'Select Audio File',...
    'E:\Users\Devon\Dropbox\02 - Projects\16 Dominos\02 - Audio Data');
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

% Crop data
xx = xx(int32(bx(1)*Fs):int32(bx(2)*Fs));
% xx = xx-xx(1);
yy = yy(int32(bx(1)*Fs):int32(bx(2)*Fs));

%% DEMON (Detection of Envelope Modulation on Noise)
% rectification
yy = sqrt(yy.^2); % Full Bridge
% yy(yy<=0) = 0;  % Half Bridge

% Envelop detection
envelopeWindow = 115;
[uppEnv, ~] = envelope(yy,envelopeWindow,'peak');

figure('Name', 'Analysed Data',...
    'OuterPosition',[0 0 screenSize(3) screenSize(4)])
subplot(2,2,1)
hold on
plot(xx,yy)
plot(xx,uppEnv,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Amp')
title('Rectified Envelope')

%% Find Peaks of Envelope
peakHeightThreshold = 0.07;
peakSeperationThreshold = 0.018;
[peakVal,peakLoc,w,prom] = findpeaks(uppEnv,Fs,...
                              'MinPeakDistance',peakSeperationThreshold,...
                              'MinPeakProminence',peakHeightThreshold);
                          
% Plot Envelop and peaks
subplot(2,2,3)
hold on
plot(xx,uppEnv)
xlabel('Time (s)')
ylabel('Amp')
title('Isolated Peaks')
hold on
plot(peakLoc(1:min(32,end))+xx(1),peakVal(1:min(32,end)),'ro')
legend('Filtered Data','Peaks')

% Time derivative of impact time (peaks)
for i=1:length(peakLoc)-1
    delLoc(i) = peakLoc(i+1)-peakLoc(i);
end

% Plot Frequency
subplot(2,2,4)
plot(peakLoc(2:min(32,end))-peakLoc(2),1./delLoc(1:min(31,end)),'b-*')
hold on
meanFreq = mean(1./delLoc(1:min(31,end)));
xl = xlim;
line(xl,meanFreq.*[1,1],'Color','r')
xlabel('Time (s)')
ylabel('Frequency (1/s)')
title(['Impact Frequency wrt Time: Mean Freq = ',num2str(meanFreq)])

%% PSD of data
[pxx,f] = periodogram(uppEnv,rectwin(length(uppEnv)),4*2^nextpow2(length(uppEnv)),Fs);
subplot(2,2,2)
hold on
plot(f,pxx)
xlabel('Frequency (1/s)')
ylabel('Amp')
title('PSD of Envelope')
xlim([5,100])

%% Save data of interest
[~,fname,~]=fileparts(fileName)
save([pathname,fname,'.mat'],'peakLoc','delLoc')
