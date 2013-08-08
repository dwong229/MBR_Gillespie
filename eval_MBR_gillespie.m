function [rtRatio,avgdtheta,dx] = eval_MBR_gillespie(timeVecAll,state,timeLim)
% use with ecoli simulation
% compute runTime, tumble Time and run/tumble time ratio
numCells = size(state.F,2);

% time is given in seconds, need to sample at the same time each round
timeVec = timeVecAll;
timeVec(timeVecAll>timeLim) = [];
% compute duration between each rxn
diffTime = timeVec(2:end) - timeVec(1:end-1);
%position
posn = state.posn(1:length(timeVec),:);
dx = ([posn(end,1:2) - posn(1,1:2)]);
theta = unwrap(deg2rad(posn(:,3)));
avgdtheta = rad2deg(theta(end)-theta(1))/timeLim;


% compute how much time each cell was in run and tumble
diffTimeMat = repmat(diffTime',[1 numCells]);
runTime = sum(state.F(1:length(diffTime),:).*diffTimeMat,1);

tumbleTime= repmat(timeVec(end),[1 numCells]) - runTime;

rtRatio = runTime./tumbleTime;

if isinf(rtRatio)
    disp('rtRatio inf')
    keyboard
end

end