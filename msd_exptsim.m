% Script for generating simulated data with experimental data updates
% User: - sets rotation or translation
%       - tau: seconds between each posn update with

clear
close all
clc 

%% USER INPUTS
transRot = 1; % 1: translation, 2: rotation
tau = 2; % amt of time to simulate between updating from expt data
frame2sec = 2; % frames per second

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
exptTime = 1/frame2sec*[0:length(th)-1];

%% Start simulation

h1 = figure('Position',[67 420 1338 537]);

subplot(1,2,1)
xysimplot = plot([0 1],[0 1],'-b');
hold on
%set(xyrealplot,'XData',xy(1:closestIdx,1),'YData',xy(1:closestIdx,2))

xyrealplot = plot(xy(:,1)*1e6,xy(:,2)*1e6,'or');
title('MBR position MSD sim with expt updates')
xlabel('X Position (um)')
ylabel('Y Position (um)')
%legend('x','y')
%xlim([0 ceil(timeVec(end))])
%axis equal

subplot(1,2,2)
orientsimplot = plot(0,0,'-b');
hold on
orientrealplot = plot(0,0,'or');
xlabel('Time (s)')
ylabel('Orientation (deg)')
axis([0,ceil(exptTime(end)),0,360])

simTime = [];
simTimeLast = 0;
mbrState = struct('posn',[],'cellposn',[],'F',[],'cellAngle',[],'detTime',[],'detPosn',[]);% simulate
while simTimeLast < exptTime(end)
    
    [~, closestIdx] = find(exptTime>simTimeLast,1,'first');
    % update mbr posn and orientation
    lastState = [xy(closestIdx,:) th(closestIdx)];
    
    % simulate for tau time
    [newTime,newmbrState] = wrapperMBRgillespiefunc(cellposn,MBRcorners,lastState,tau);
    
    newTimeLength = find(newTime<tau,1,'last');
    
    simTime = cat(2,simTime,newTime(1:newTimeLength)+exptTime(closestIdx));
    mbrState.posn = cat(1,mbrState.posn,newmbrState.posn(1:newTimeLength,:));
    
    simTimeLast = simTime(end);
    
    %% plot update
    %plot simulation
    x = mbrState.posn(:,1);
    y = mbrState.posn(:,2);
    thSim = mbrState.posn(:,3);
    
    %subplot(2,1,1)
    
    set(xysimplot,'XData',x,'YData',y)
    %set(xyrealplot,'XData',xy(1:closestIdx,1),'YData',xy(1:closestIdx,2))
    set(orientsimplot,'XData',simTime,'YData',thSim)
    set(orientrealplot,'XData',exptTime(1:closestIdx),'YData',th(1:closestIdx))
    
    keyboard
end

