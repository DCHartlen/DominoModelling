%% Shaw Domino Model
% 
% Implimented based on van Leeuwen (2004). The Domino Effect, ArXiv Physics
% E-Prints.
%
% Created by:  D.C. Hartlen, EIT
% Date:        06-Apr-2018
% Modified by: D.C. Hartlen, EIT 
% Date:        07-May-2018
%
% Change log:
% 07-May-2018
% Adapted code to allow for arbitrary domino spacing. Previous code
% enforced constant domino spacing, preventing the initial domino spacing
% to be increased to allow for consistant initiation energy.

%% Initialization
close all
clear
clc

%% Domino and Physical Parameters
NDom = 50;   % Number of dominoes in chain
mDom = 8e-3; % Mass of dominoes
% a = 8e-3;   % Width of a single domino
% b = 50e-3;  % Height of a single domino
% c = 30e-3;  % Chain spacing. This includes the width of one domino
a = 7.45e-3;   % Width of a single domino
b = 48e-3;  % Height of a single domino
c = (1+2)*a;  % Chain spacing. This includes the width of one domino
cInit = (1+3)*a; % Spacing of first domino. Used to get more energy in system.
I = mDom*(a^2+b^2)/3; % Domino's mass moment of inertia about corner
mu = 0.15;  % Coefficient of friction

g = 9.81;   % Acceleration of gravity
initAngle = deg2rad(8.9); % Angle which first domino starts at
minTipAngle = asin(a/2/sqrt((0.5*a)^2+(0.5*b)^2));

%% Time integration parameters
nASteps = 50;

%% Chain Initialization
theta = 0.*ones(NDom,1); % All dominoes start vertical
theta(1) = initAngle; % Start first domino at 20deg from vertical 
omega = 1.0e-10.*ones(NDom,1); %Initialize speed with very small number (avoids NaN Issues)
time = 0.*ones(NDom,1); % Initalize fall time for each domino
cArray = c.*ones(NDom,1); % Initialize array of spacings
cArray(1) = cInit;  % Set first element of spacing array to be initial spacin

thetaCrit = asin((c-a)/b); % Critical angle at contact
thetaCritInit = asin((cInit-a)/b); % Critical angle at contact for first domino

%% Initiate integration level output files
nTimesteps = NDom*nASteps;
thetaOut = zeros(nTimesteps,NDom); 
omegaOut = zeros(nTimesteps,NDom);
energyTotOut = zeros(nTimesteps,1);
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
    energyTotOut((nDom-1)*nASteps+i) = Ei;
    timeOut((nDom-1)*nASteps+i) = tn;
end
theta(nDom) = theta(nDom);
omega(nDom) = omega(nDom);
time(nDom) = tn;

%% Start loop for remaining dominoes

for nDom = 2:NDom
%     nDom = 2
    fprintf('--------------------Domino %d--------------------\n',nDom)
    % Start van Leewen contact modeling (accounts for changing moment arms)
    % Compute alpha and beta coefficients as well as r. This calculuation
    % is done back to front, instead of the reversed arrays documented by
    % van Leewen.
    r = zeros(nDom,1);
    r(end) = 1;
    for iRec = nDom-1:-1:1
        alpha(iRec) = b*cos(theta(iRec)-theta(iRec+1)) -...
            cArray(iRec)*sin(theta(iRec+1)) - mu*a;
        beta(iRec) = b*cos(theta(iRec)-theta(iRec+1)) + ...
            mu*b*sin(theta(iRec)-theta(iRec+1));
        r(iRec) = r(iRec+1)*alpha(iRec)/beta(iRec);
    end
    W = [W;1];
    J = sum(r.*W);
    omega(nDom) = (J-1)/J*omega(nDom-1);

    % Begin recusion to compute speed of previous dominoes after impact
    for iRec = nDom-1:-1:1
        omega(iRec) = omega(iRec+1)*...
            (1- (cArray(iRec)*sin(theta(iRec+1)))/...
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
                asin(cArray(iRec)/b*cos(theta(iRec+1))-a/b);
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
                (1- ((cArray(iRec)*sin(theta(iRec+1)))/...
                (b*cos(theta(iRec)-theta(iRec+1)))));
        end
        
        % Iteration energy balance for a check
        Ei = 0.5*mDom*g*(b*sum(cos(theta(1:nDom)))+a*sum(sin(theta(1:nDom))))...
            + 0.5*I*sum(omega(1:nDom).^2);
        
        % Write to output matrices
        thetaOut((nDom-1)*nASteps+i,:) = theta';
        omegaOut((nDom-1)*nASteps+i,:) = omega';
        energyTotOut((nDom-1)*nASteps+i) = Ei;
        timeOut((nDom-1)*nASteps+i) = sum(time(1:nDom-1))+tn;
         
    end % End angle integration
    
    % Store variables from this domino. This is not a time history
    theta(nDom) = theta(nDom);
    omega(nDom) = omega(nDom);
    time(nDom) = tn;
end
velocity = cArray(iRec)./time;
velocityND = velocity./sqrt(g*b);