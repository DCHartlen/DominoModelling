close all
clear
clc

%% Crop data
screenSize = get( groot, 'Screensize' );

% Load Data
[fileNames, pathname, ~] = uigetfile({'*.*'},'Select Audio File',...
    'E:\Users\Devon\Dropbox\02 - Projects\16 Dominos\02 - Audio Data',...
    'MultiSelect', 'on');

for i = 1:length(fileNames)
    load([pathname,fileNames{i}]);
    Data(i).peakLoc = peakLoc;
    Data(i).delLoc = delLoc;
end

figure()
hold on
for i = 1:length(fileNames)
    plot(Data(i).peakLoc(2:min(32,end))-Data(i).peakLoc(2),...
        1./Data(i).delLoc(1:min(31,end)),'-*')
end
legend('1','2','3','4','5')
xlabel('Time of Impact (s)')
ylabel('ImpactFreq (1/s)')

figure()
hold on
for i = [1,2,3,4,5]
    plot([1:1:31],...
        1./Data(i).delLoc(1:min(31,end)),'-*')
end
legend('1','2','3','4','5')
xlabel('Domino Number')
ylabel('ImpactFreq (1/s)')

for i=1:length(Data)
    freqData(i,:) = 1./Data(i).delLoc(1:min(31,end));
end

tscore = tinv(0.95,length(Data)-1)
for i=1:length(freqData)
    meanFreq(i) = mean(freqData(:,i));
    stdFreq(i) = std(freqData(:,i));
        
    uncertainty(i) = tscore*stdFreq(i)/sqrt(length(Data));
end

figure()
hold on
errorbar([1:1:31],meanFreq,uncertainty)
% xlabel('Domino Number')
% ylabel('ImpactFreq (1/s)')
% for i = [1,2,3,4,5]
%     plot([1:1:31],...
%         1./Data(i).delLoc(1:min(31,end)),'-*')
% 
% end
% legend('1','2','3','4','5')




