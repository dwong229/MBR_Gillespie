%% Compute MSD for tracked MBRs
close all
clear all
load('2Htrackstat.mat')
% H2angle (in deg)
% H2xy (in microns)

x = H2xy(:,1);
y = H2xy(:,2);

tau = 1;

figure('Position',[114 461 1230 454]);
subplot(1,2,1)
for tau = 1:20
    msd{tau} = sum((H2xy(end:-1:1+tau,:) - H2xy(end-tau:-1:1,:)).^2,2);
    meanval(tau) = mean(msd{tau});
    stddev(tau) = std(msd{tau});

    hold on
    plot(tau,msd{tau},'ob')
    
end

plot(1:tau,meanval,'s-r')
errorbar(1:tau,meanval,stddev,'r')


xlabel('\tau')
ylabel('msd(\tau)')
title('Translating: H2 dataset')

H2 = struct('xy',H2xy,'msd',msd,'meanval',meanval,'stddev',stddev);


load('H3trackstat.mat')
% H3angle (in deg)
% H3xy (in microns)
subplot(1,2,2)
for tau = 1:20
    msd{tau} = sum((H3xy(end:-1:1+tau,:) - H3xy(end-tau:-1:1,:)).^2,2);
    meanval(tau) = mean(msd{tau});
    stddev(tau) = std(msd{tau});
    hold on
    plot(tau,msd{tau},'ob')
    
end

plot(1:tau,meanval,'s-r')
errorbar(1:tau,meanval,stddev,'r')

xlabel('\tau')
ylabel('msd(\tau)')
title('Rotating: H3 dataset')
H3 = struct('xy',H3xy,'msd',msd,'meanval',meanval,'stddev',stddev);

print(gcf,'-djpeg','msd.jpg')


%% 
figure('Position',[114 461 1230 454])
ax(1) = subplot(1,2,1);
%hold on
loglog(1:tau,H2(1).meanval,'s-r')
%keyboard
hold on
for i = 1:tau
    loglog(i,H2(i).msd,'ob');
end


xlabel('\tau')
ylabel('msd(\tau)')
title('Translating: H2 dataset')

ax(2) = subplot(1,2,2);
loglog(1:tau,H3(1).meanval,'s-r')
hold on
for i = 1:tau
    loglog(i,H3(i).msd,'ob');
end

xlabel('\tau')
ylabel('msd(\tau)')
title('Rotating: H3 dataset')
linkaxes([ax(1),ax(2)],'xy')
print(gcf,'-djpeg','msdloglog.jpg')