disp('Running H3 rotating comparison plots')

load('H3trackstat.mat')
% H3angle (in deg)
% H3xy (in microns)


% determine when to stop simualtion
endIdx = find(timeVec>12.5,1);
if isempty(endIdx)
    endIdx = length(timeVec);
end

% plot trajectory
htraj = figure
plot(H3xy(:,1),H3xy(:,2),'.r')
hold on
plot(state.posn(1:endIdx,1),state.posn(1:endIdx,2),'-g')
legend('Experiment','Simulation')
xlabel('x (um)')
ylabel('y (um)')

axis equal
axis ij

htheta = figure
timeVecTrack = 1/4*[0:length(H3angle)-1];
H3angle = -(H3angle); %flip into image coords

plot(timeVecTrack,H3angle,'.r')
hold on
state.posn(:,3) = rad2deg(unwrap(deg2rad(state.posn(:,3)))) - state.posn(1,3);
state.posn(:,3) = [state.posn(:,3)>360]*(-360) + state.posn(:,3);
plot(timeVec(1:endIdx),state.posn(1:endIdx,3),'-g')
legend('Experiment','Simulation')
axis([0,max(timeVecTrack(end),timeVec(end)),0,360])
ylabel('Orientation (deg)')
xlabel('Time (sec')