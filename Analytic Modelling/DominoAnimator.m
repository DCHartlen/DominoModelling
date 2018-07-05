function DominoAnimator(dominoDims, spacing, timestamps, dominoAngles)
%% Domino Plotter
% 
% Takes angle and postion data and plots dominos falling
%
% Created by:  P.C. Hennessey, EIT
% Date:        09-Apr-2018
% Modified by: D.C. Hartlen, EIT
% Date:        05-Jul-2018
%
% INPUTS:
%   dominoDims: array of the format [height, width] of domino.
%   spacing: array containing all arbitrary spacing
%   timestamp: array of timestamps corresponding to dominoAngles
%   dominoAngles: contains angle data for all dominoes over all timesteps.
%
% TODO: Add ability to cope with arbitrary spacing
% TODO: Impliment better data input method which does not rely on 
%       variables being in the workspace already. Turn into a function?

close all

%% Assign parameters 

% extract domino dimensions for ease of use
w = dominoDims(2)*1000; %Domino Width
h = dominoDims(1)*1000; %Domino Height

s = 4*w; %Domino Spacing

% extract number of dominoes
nDom = size(dominoAngles,2); %number of dominos
nTimestamp = length(timestamps); %number timesteps

tmax = max(timestamps);

%% Create images for the animation
% Outer loop creates snapshot of domino position
for  i = 1:10:nTimestamp
    % Inner loop creates the dominos one by one
    for j = 1:nDom
        % TODO Update spacing logic
        d=j*s;
        theta = pi/2 - dominoAngles(i,j);
        %Setup point equations
        pt1 = [d; 0];
        pt2 = [d+cos(theta)*h; sin(theta)*h];
        pt3 = [d+cos(theta)*h-sin(theta)*w; sin(theta)*h+cos(theta)*w];
        pt4 = [d-sin(theta)*w; cos(theta)*w];
        
        domPts = [pt1, pt2, pt3, pt4];
        
        % TODO Find a more elegant plotting? Maybe call figure outside
        % inner loop?
        figure(1)
        title('Shaw Model (Analytic) 25 Domino Numeric Simulation')
        set(gcf, 'Position', [100, 100, 1000, 100])
        fill(domPts(1,:), domPts(2,:), 'r')
        xlabel('Length (mm)')
        ylabel('Height (mm)')
        axis([-5,s*nDom,-5,h+5])
        drawnow limitrate
        hold on
    end
    
    % Convert a figure to an image and store
    frame = getframe(figure(1));
    im{i} = frame2im(frame);
    [A,map] = rgb2ind(im{i},256);
        
    hold off
end

%% Write the images to file
filename = 'testGIF.gif';

for idx = 1:50:nTimestamp
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.002);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.002);
    end
end




