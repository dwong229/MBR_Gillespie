%% Compute MSD for tracked MBRs
close all
clear all
load('2Htrackstat.mat')
% H2angle (in deg)
% H2xy (in microns)

x = H2xy(:,1);
y = H2xy(:,2);

sqrtmsd = false;
sqrtt = false;
tRange = 1:20;

%% MSD-time plot
hfigmsd(1) = figure('Position',[114 461 1230 454]);
% H2
H2 = computemsd(tRange,H2xy,sqrtmsd,sqrtt);

subplot(1,2,1);
% plot trials
for idx = 1:length(tRange)
    hold on
    plot(H2.tRange(idx),H2.msd{idx},'ob')
end

% means values with std
hfigmsd(1) = plot(H2.tRange,H2.meanval,'s-r');
errorbar(H2.tRange,H2.meanval,H2.stddev,'r','LineWidth',5)

xlabel('\tau')
if sqrtmsd
    xlabel('sqrt(\tau)')
end
ylabel('msd(\tau)')
if sqrtmsd
    ylabel('sqrt(msd(\tau))')
end
title('Translating: H2 dataset')

% H3
load('H3trackstat.mat')
% H3angle (in deg)
% H3xy (in microns)
hfigmsd(2) = subplot(1,2,2)
H3 = computemsd(tRange,H3xy,sqrtmsd,sqrtt);

% plot trials
for idx = 1:length(tRange)
    hold on
    plot(H3.tRange(idx),H3.msd{idx},'ob')
end

% means values with std
plot(H3.tRange,H3.meanval,'s-r')
errorbar(H3.tRange,H3.meanval,H3.stddev,'r','LineWidth',5)

plot(H3.tRange,H3.meanval,'s-r')
errorbar(H3.tRange,H3.meanval,H3.stddev,'r','LineWidth',5)

xlabel('\tau')
if sqrtmsd
    xlabel('sqrt(\tau)')
end
ylabel('msd(\tau)')
if sqrtmsd
    ylabel('sqrt(msd(\tau))')
end
title('Rotating: H3 dataset')

%print(gcf,'-djpeg','msd.jpg')

%% Diffusion plot
%figure(hfigmsd(1))
subplot(1,2,1)
q = 4; % numericalconstant which depends on dimensionality.
D = H2.meanval(1);
ydiff = D*[0 tRange];
plot([0 tRange],ydiff,'-g','LineWidth',5)

subplot(1,2,2)
q = 4; % numericalconstant which depends on dimensionality.
D = H3.meanval(1);
ydiff = D*[0 tRange];
plot([0 tRange],ydiff,'-g','LineWidth',5)

%print(gcf,'-djpeg','msd_diffusioncompare.jpg')

keyboard

%% Log-Log Graph
figure('Position',[114 461 1230 454])
ax(1) = subplot(1,2,1);
%hold on
loglog(H2.tRange,H2.meanval,'s-r')
%keyboard
hold on
for i = 1:length(tRange)
    loglog(H2.tRange(i),H2.msd{i},'ob');
end


xlabel('\tau')
if sqrtmsd
    xlabel('sqrt(\tau)')
end
ylabel('msd(\tau)')
if sqrtmsd
    ylabel('sqrt(msd(\tau))')
end
title('LOGLOG Translating: H2 dataset')

ax(2) = subplot(1,2,2);
loglog(H3.tRange,H3.meanval,'s-r')
hold on
for i = 1:length(tRange)
    loglog(H3.tRange(i),H3.msd{i},'ob');
end

xlabel('\tau')
if sqrtmsd
    xlabel('sqrt(\tau)')
end
ylabel('msd(\tau)')
if sqrtmsd
    ylabel('sqrt(msd(\tau))')
end
title('LOGLOG Rotating: H3 dataset')
linkaxes([ax(1),ax(2)],'xy')
%print(gcf,'-djpeg','msdloglog.jpg')