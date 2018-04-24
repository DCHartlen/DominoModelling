%% Shaw Domino Model
% 
% Implimented based on Shaw (1978). Mechanics of a chain of dominoes.
% American Journal of Physics 46(6). pp 640-642.
%
% Created by:  D.C. Hartlen, EIT
% Date:        04-Apr-2018
% Modified by:  
% Date:        

%% Initialization
close all
clear
clc

%% Domino and Physical Parameters
NDom = 32;   % Number of dominoes in chain
mDom = 8e-3; % Mass of dominoes
% a = 8e-3;   % Width of a single domino
% b = 50e-3;  % Height of a single domino
% c = 30e-3;  % Chain spacing. This includes the width of one domino
a = 7.45e-3;   % Width of a single domino
b = 48e-3;  % Height of a single domino
c = (1+2)*a;  % Chain spacing. This includes the width of one domino
cInit = (1+3)*a; % Spacing of first domino. Used to get more energy in system.
I = mDom*(a^2+b^2)/3; % Domino's mass moment of inertia about corner

g = 9.81;   % Acceleration of gravity
initAngle = deg2rad(10); % Angle which first domino starts at
minTipAngle = asin(a/2/sqrt((0.5*a)^2+(0.5*b)^2));

%% Time integration parameters
nASteps = 50;

%% Chain Initialization
theta = 0.*ones(NDom,1); % All dominoes start vertical
theta(1) = initAngle; % Start first domino at 20deg from vertical 
omega = 1.0e-10.*ones(NDom,1); %Initialize speed with very small number (avoids NaN Issues)
time = 0.*ones(NDom,1); % Initalize fall time for each domino

thetaCrit = asin((c-a)/b); % Critical angle at contact
thetaCritInit = asin((cInit-a)/b); % Critical angle at contact for first domino

%% Initiate integration level output files
nTimesteps = NDom*nASteps;
thetaOut = zeros(nTimesteps,NDom); 
omegaOut = zeros(nTimesteps,NDom);
energyOut = zeros(nTimesteps,1);
timeOut = zeros(nTimesteps,1);

%% Error check of input parameters
% Ensure initiation angle is great enough to self-start (CoG outside
% domino)
assert(theta(1)>minTipAngle,...
    'Error: Initiation angle does not provide self-tip')
% Ensure initiation angle does not exceed critical angle
assert(theta(1)<thetaCritInit, 'Error: Initiation angle is larger than critical angle %d',rad2deg(thetaCritInit))

%% Start Calculation for domino 1
% Pre-integration compuation
nDom = 1;  % Special case of at first domino
fprintf('--------------------Domino %d--------------------\n',nDom)

thetaInit = theta(nDom); % Initial angle of domino n
delTheta = (thetaCritInit-thetaInit)/nASteps; % Compute angle step for integration
% Energy in first n dominoes (energy in first domino)
E = 0.5*mDom*g*(b*sum(cos(theta(1:nDom)))+a*sum(sin(theta(1:nDom))))...
    + 0.5*I*sum(omega(1:nDom).^2)
tn = 0; %time counter for current domino

% Commence Angle/Time Integration
for i=1:nASteps
    % Add angle step to domino N
    theta(nDom) = theta(nDom)+delTheta;
    
    % set up recursive coefficients
    U = cos(theta(1:nDom))./cos(theta(nDom));
    V = sin(theta(1:nDom))./sin(theta(nDom));
    W = omega(1:nDom)./omega(nDom);
    
    % Compute angular velocity of Domino N from energy and angle
    omega(nDom) = ((2*E-mDom*g*(b*cos(theta(nDom))*sum(U)+a*sin(theta(nDom))*sum(V))) / ...
        (I*sum(W.^2)))^0.5;
    % Use incremental angle and instantaneous velocity to compute time
    % increment
    tn = tn + delTheta/omega(nDom);
    
    % Iteration energy balance for a check
    Ei = 0.5*mDom*g*(b*sum(cos(theta(1:nDom)))+a*sum(sin(theta(1:nDom))))...
    + 0.5*I*sum(omega(1:nDom).^2);
    
    % Write to output matrices
    thetaOut((nDom-1)*nASteps+i,:) = theta';
    omegaOut((nDom-1)*nASteps+i,:) = omega';
    energyOut((nDom-1)*nASteps+i) = Ei;
    timeOut((nDom-1)*nASteps+i) = tn;
end
theta(nDom) = theta(nDom);
omega(nDom) = omega(nDom);
time(nDom) = tn;

%% Start loop for remaining dominoes

for nDom = 2:NDom
    fprintf('--------------------Domino %d--------------------\n',nDom)
    % Impact rule for dominoes
    % Find speed of next domino in line
    omega(nDom) = omega(nDom-1)*(sum(W)/sum([W;1]));
    % Begin recusion to compute speed of previous dominoes after impact
    for iRec = nDom-1:-1:1
        omega(iRec) = omega(iRec+1)*...
            (1- (c*sin(theta(iRec+1)))/...
            (b*cos(theta(iRec)-theta(iRec+1))));
    end
    
    % Calculate active energy
    E = 0.5*mDom*g*(b*sum(cos(theta(1:nDom)))+a*sum(sin(theta(1:nDom))))...
        + 0.5*I*sum(omega(1:nDom).^2)
    
    % Start Angle/Time integration
    thetaInit = theta(nDom); % Initial angle of domino n
    delTheta = (thetaCrit-thetaInit)/nASteps; % Compute angle step for integration
    tn = 0; %time counter for current domino
    
    for i = 1:nASteps
        % Add angle step to domino N
        theta(nDom) = theta(nDom)+delTheta;
        
        % Run recursion to compute angle of previous dominoes after
        % incremental step
        for iRec = nDom-1:-1:1
            theta(iRec) = theta(iRec+1) +...
                asin(c/b*cos(theta(iRec+1))-a/b);
        end
        
        % set up recursive coefficients
        U = cos(theta(1:nDom))./cos(theta(nDom));
        V = sin(theta(1:nDom))./sin(theta(nDom));
        W = omega(1:nDom)./omega(nDom);
        
        % Compute angular velocity of Domino N from energy and angle
        omega(nDom) = ((2*E-mDom*g*(b*cos(theta(nDom))*sum(U)+a*sin(theta(nDom))*sum(V))) / ...
            (I*sum(W.^2)))^0.5;
        assert(isreal(omega(nDom)),'Error: insufficient energy to topple domino %i',nDom)
        % Use incremental angle and instantaneous velocity to compute time
        % increment
        tn = tn + delTheta/omega(nDom);
        
        % Run recusion to compute speed of all previous dominoes
        for iRec = nDom-1:-1:1
            omega(iRec) = omega(iRec+1)*...
                (1- ((c*sin(theta(iRec+1)))/...
                (b*cos(theta(iRec)-theta(iRec+1)))));
        end
        
        % Iteration energy balance for a check
        Ei = 0.5*mDom*g*(b*sum(cos(theta(1:nDom)))+a*sum(sin(theta(1:nDom))))...
            + 0.5*I*sum(omega(1:nDom).^2);
        
        % Write to output matrices
        thetaOut((nDom-1)*nASteps+i,:) = theta';
        omegaOut((nDom-1)*nASteps+i,:) = omega';
        energyOut((nDom-1)*nASteps+i) = Ei;
        timeOut((nDom-1)*nASteps+i) = sum(time(1:nDom-1))+tn;
         
    end % End angle integration
    
    % Store variables from this domino. This is not a time history
    theta(nDom) = theta(nDom);
    omega(nDom) = omega(nDom);
    time(nDom) = tn;
end

velocity = c./time;
velocityND = velocity./sqrt(g*b);
