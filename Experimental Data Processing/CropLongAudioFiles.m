%% Crop Long Recordings
%
% The process audio files works best with sound files which consist
% primarily of the domino noises. This script allows a user to
% interactively crop a long audio file and save the cropped data as a
% smaller one. 
%
% Created by:  D.C. Hartlen, EIT
% Date:        17-Apr-2018
% Modified by: D.C. Hartlen, EIT
% Date:        17-Aug-2018

close all
clear
clc

screenSize = get( groot, 'Screensize' );

% Load Data
% [yy,Fs] = audioread('Domino_7pt5_1.m4a');
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

[path, name, ext] = fileparts(fileName);
name = [name '_Cropped'];
newFileName = [name,'.wav'];

audiowrite([pathname newFileName], yy,Fs)