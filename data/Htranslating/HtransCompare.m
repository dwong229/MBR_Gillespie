disp('Running H translating comparison plots')

load('2Htrackstat.mat')
% H2angle (in deg)
% H2xy (in microns)


% determine when to stop simualtion
endIdx = find(timeVec>7,1);
if isempty(endIdx)
    endIdx = length(timeVec);
end

% plot trajectory
htraj = figure
title('Translating H: Comparison to experimental data')
plot(H2xy(:,1),H2xy(:,2),'.r')
hold on
subplot(1,2,2)
plot(state.posn(:,1),state.posn(:,2),'-g')
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
state.posn(:,3) = state.posn(:,3) - state.posn(1,3);
plot(timeVec,state.posn(:,3),'-g')
legend('Experiment','Simulation')
ylabel('Orientation (deg)')
xlabel('Time (sec')

% Comparison to non-stochastic results

