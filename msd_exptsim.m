% Script for generating simulated data with experimental data updates
% User: - sets rotation or translation
%       - tau: seconds between each posn update with

%% USER INPUTS
transRot = 1; % 1: translation, 2: rotation
tau = 2; % amt of time to simulate between updating from expt data

%% Load appropriate files for translation or rotation mode
if transRot == 1
    disp('Translation')
    load('2Htrackstat.mat')
    xy = H2xy; % nx2
    th = H2angle; % in deg
    
    % load cell posn:
    cellposnfile = 'cellposnOpenCV2H_headangle.mat';
    
    % load MBRcorners
    % H 60 x 60
    % translating H
    MBRcorners.cells(:,1) = [-30;25]; %x coordinates
    MBRcorners.cells(:,2) = [-30;25]; %y coordinates
    
    MBRcorners.nocells = [-18 8;13 30;...
        -12 -30;13 -5];
else
    disp('Rotation')
    load('H3taackstat.mat')
    xy = H3xy; % nx2
    th = H3angle; % in deg
    
    % load cell posn:
    cellposnfile = 'headangle_data_H3reverse.mat';
    
    % load MBRcorners
    % H 60 x 60
    % rotating H
    MBRcorners.cells(:,1) = [-30;30]; %x coordinates
    MBRcorners.cells(:,2) = [-30;30]; %y coordinates
    
    MBRcorners.nocells = [-18 8;18 30;-18 -30;18 -8];
end

load(cellposnfile)

%% Start simulation

frame2sec = 2; % frames per second

exptTime = frame2sec*[0:length(th)-1];

simTimeLast = 0;
% simulate
while simTime < exptTime(end)
    
    [~, closestIdx] = min(abs(exptTime - simTime));
    % update mbr posn and orientation
    lastState = [xy(closestIdx,:) th(closestIdx)];
    
    % simulate for tau time
    [newTime,newmbrState] = wrapperMBRgillespiefunc(cellposn,MBRcorners,lastState,tau);
    
    [~,newTimeLength] = min(abs(newTime - tau));
    
    simTime = cat(simTime,newTime(1:newTimeLength)+simTimeLast);
    mbrState = cat(1,mbrState,newmbrState);
    
    
    simTimeLast = simTime(end);
    
    % plot update
    keyboard
end