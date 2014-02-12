function [] = HtransCompare(MBRx,MBRy,MBRth,timeVec,varargin)
 
disp('Running H translating comparison plots')

if nargin == 4
    plotstr = '-y';
elseif varargin{1} == 'det'
    plotstr = '-b';
elseif varargin{1} == 'sim'
    plotstr = '-g';
end

load('2Htrackstat.mat')
% H2angle (in deg)
% H2xy (in microns)


% determine when to stop simualtion
endIdx = find(timeVec>7,1);
endIdx = length(timeVec);
if isempty(endIdx)
    endIdx = length(timeVec);
end

% plot trajectory
htraj = figure
title('Translating H: Comparison to experimental data')
subplot(1,2,1)
plot(H2xy(:,1),H2xy(:,2),'.r')
hold on
plot(MBRx,MBRy,plotstr)
legend('Experiment','Simulation')
xlabel('x (um)')
ylabel('y (um)')
axis equal
axis ij

title('Translating H: Comparison to experimental data')
timeVecTrack = 1/5*[0:length(H2angle)-1];
H3angle = -(H2angle - H2angle(1)); %flip into image coords
subplot(1,2,2)
plot(timeVecTrack,H3angle,'.r')
hold on
MBRth = MBRth - MBRth(1);
plot(timeVec,MBRth,plotstr)
legend('Experiment','Simulation')
ylabel('Orientation (deg)')
xlabel('Time (sec')

% Comparison to non-stochastic results

