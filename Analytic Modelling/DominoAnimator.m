%% Domino Plotter
% 
% Takes angle and postion data and plots dominos falling
%
% Created by:  P.C. Hennessey, EIT
% Date:        09-Apr-2018
% Modified by: 
% Date:        
%
% TODO: Add ability to cope with arbitrary spacing
% TODO: Impliment better data input method which does not rely on 
%       variables being in the workspace already. Turn into a function?
close all


%PROGRAM INPUTS
%Define domino size and spacing (mm)
%load domino data 
%input thetaOut and timeOut

w = 7.45; %Domino Width
h = 48; %Domino Height
s = 4*w; %Domino Spacing
n = 32; %number of dominos
t = 1600; %number timesteps

tmax = max(timeOut);

%Build Matrix of plot pts

for  i = 1:10:t

    for j = 1:n
        
        d=j*s;
        theta = pi/2 - thetaOut(i,j);
        %Setup point equations
        pt1 = [d; 0];
        pt2 = [d+cos(theta)*h; sin(theta)*h];
        pt3 = [d+cos(theta)*h-sin(theta)*w; sin(theta)*h+cos(theta)*w];
        pt4 = [d-sin(theta)*w; cos(theta)*w];
        
        domPts = [pt1, pt2, pt3, pt4];
        
        figure(1)
        xlabel('Length (mm)')
        ylabel('Height (mm)')
        title('Shaw Model (Analytic) 25 Domino Numeric Simulation')
        set(gcf, 'Position', [100, 100, 1000, 100])
        fill(domPts(1,:), domPts(2,:), 'r')
        axis([-5,s*n,-5,h+5])
        drawnow limitrate
        hold on
    end
    
    frame = getframe(figure(1));
    im{i} = frame2im(frame);
    [A,map] = rgb2ind(im{i},256);
        
    hold off
end
%%
filename = 'testGIF.gif';

for idx = 1:50:t
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.002);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.002);
    end
end




