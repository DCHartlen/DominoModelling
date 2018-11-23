close all
clear
clc

%% Crop data
screenSize = get( groot, 'Screensize' );

% Load Data
[fileNames, pathname, ~] = uigetfile({'*.*'},'Select Audio File',...
    'E:\Users\Devon\Dropbox\02 - Projects\16 Dominos\02 - Audio Data',...
    'MultiSelect', 'on');

figure()
hold on
for i=1:length(fileNames)
    load([pathname,fileNames{i}]);
    errorbar([1:1:31],meanFreq,uncertainty)
end
xlabel('Domino Number')
ylabel('Impact Frequency')
legend('1 Domino Spacing','2 Domino Spacing','3 Domino Spacing')