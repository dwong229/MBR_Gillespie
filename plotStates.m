function [] = plotStates(timeVec,state)

% Plot trajectory
%plotstate = zeros(size(state));
for i =1:size(state.chem,2)
    %[plottime,plotstate(:,i)] = stairs(timeVec,state.chem(:,i));
    [plottime,plotstate(:,i)] = stairs(timeVec,state.chem(:,i));
end

chem = figure
subplot(4,1,1) %I,A
plot(plottime,plotstate(:,1),'r',plottime,plotstate(:,2),'g')
title('Inactive/Active')

subplot(4,1,2) %CheA and CheAp
plot(plottime,plotstate(:,3),'b',plottime,plotstate(:,4),'r')
%plot(plottime,plotstate(:,4),'r')
title('AA/AAp')

subplot(4,1,3) %CheY and CheYp
plot(plottime,plotstate(:,5),'b',plottime,plotstate(:,6),'r')
%plot(plottime,plotstate(:,6),'r')
title('YY/YYp')

subplot(4,1,4) %Mot
plot(plottime,plotstate(:,7),'b',plottime,plotstate(:,8),'-r')
%plot(plottime,plotstate(:,8),'-r','LineWidth',2)
title('Mot/Run')

%% 
disp = figure
displacement = [state.dyn(2:end,1:2)-state.dyn(1:end-1,1:2)];
mag = sqrt(sum(displacement.^2,2));
ang = state.dyn(1:end-1,3);
bin = 0:90:360;
angMag = zeros(size(bin));
for i = 1:numel(ang)
    binNo = find(ang(i)>bin,1,'last');
    angMag(binNo) = angMag(binNo) + mag(i);
end
normangmag = angMag/sum(angMag);
binc = bin(1:end-1)+45;
bar(binc,normangmag(1:end-1))


%% 
hist = figure
plot(state.dyn(:,1),state.dyn(:,2),'-k','LineWidth',3)
hold on
plot(state.dyn(:,1),state.dyn(:,2),'.b','MarkerSize',5)
plot(state.dyn(1,1),state.dyn(1,2),'.g','MarkerSize',20)
plot(state.dyn(end,1),state.dyn(end,2),'.r','MarkerSize',20)
xlabel('x')
ylabel('y')
title('E.coli position')

%% save plots
nowstr = datestr(now,13);
nowstr(nowstr==':') = [];
%save(cat('disp',nowstr)
