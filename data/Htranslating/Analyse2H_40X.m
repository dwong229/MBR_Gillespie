% analyse 2H_40X.avi from file full of corners
clear
close all
fps = 10;

% extract file of corners in pixels
load('FinalCornersALL2H_40X.mat')

nFrames = 200;

for i = 1:170
    idx = (i-1)*4 + 1;
    
    % find mean by computing centroid of 4 corners
    % [xmean ymean]
    centroid(i,:) = mean(corners(idx:idx+3,:),1);
    
    % find theta by computing rotation matrix between corners of current
    % and last frame
    
    %given inliers, compute r and t
    
    if i > 1
        
        % put current 
        pNow0 = [corners(idx:idx+3,:) ones(4,1)]';
        
        uNow = mean(pNow0,2);
        
        pLast = bsxfun(@minus,pLast0,uLast);
        pNow = bsxfun(@minus,pNow0,uNow);
        
        G = pNow*pLast';
        [U D V] = svd(G);
        
        S = eye(3);
        if det(G)<0
            S(3,3) = -1;
        end
        
        rInliers= U*S*V';
        tInliers= uLast - rInliers*uNow;
                
        [roll,pitch,yaw] = RotToRPY_ZXY(rInliers); %in image frame
        dtheta(i) = rad2deg(yaw);
    end
    pLast0 = [corners(idx:idx+3,:) ones(4,1)]';
    uLast = [centroid(i,:) 1]';
    
    for j = 1:4
        % organize corners so that they correspond
        % take current corner and set to last
        lastcorner(j,:) = corners(idx+j-1,:);
        
        % find next corner that is closest to lastcorner
        distMat = dist([lastcorner(j,:);corners(idx+4:idx+7,:)]');
        [a minIdx] = min(distMat(1,2:end));
        cornersTemp(j,:) = corners(idx+3+minIdx,:);
        
    end
    corners(idx+4:idx+7,:) = cornersTemp;
    lastcorner = [];
    cornersTemp = [];
end

timeVec = [1:i]*1/fps;
subplot(2,1,1)
theta = cumsum(dtheta);
plot(timeVec,theta,'.b')
title('World Frame Angle')
xlabel('Time (sec)')
ylabel('Angle (deg)')

% convert pixels to microns
% [267,760] to [283 316]
pix2mic = 494/60; % pixels per 60um

centroidMicron = bsxfun(@minus,centroid,centroid(1,:))/pix2mic;
subplot(2,1,2)
title('World Frame Position')
hold on
plot(timeVec,centroidMicron(:,1),'*b');
plot(timeVec,centroidMicron(:,2),'*g');
legend('X','Y')
xlabel('Time (sec)')
ylabel('Position (deg)')

%% rotote translation to Body Frame
% lots of cells are aligned along body frame 
pts(1,:) = [346 740]; % bottom left H
pts(2,:) = [269 345];% top left H

initialBodyAngle = rad2deg(atan2(abs(diff(pts(:,2))),diff(pts(:,1)))) - 90;

% Add initial body angle to (Angle to bottom of H from horizontal)
bodyTheta = theta + initialBodyAngle;

dxdy= diff(centroidMicron);

for k = 1:length(dxdy)
    R = @(x)[cosd(x) -sind(x);sind(x) cosd(x)]; % Body = R*World
    dxdyBody(k,:) = R(bodyTheta(k))*dxdy(k,:)';
end

dxdyBody = dxdyBody;

figure
title('Body Frame Velocity')
hold on
plot(timeVec(2:end),dxdyBody(:,1),'-b');
plot(timeVec(2:end),dxdyBody(:,2),'-g');
legend('X','Y')
xlabel('Time (sec)')
ylabel('Velocity (um/s)')

load('2H_40x.mat')
mbrStat = struct('centroid',zeros(1,2,'double'),'orientation',zeros(1,'double'),'corners',zeros(4,2,'double'));
for i = 1:170
    mbrStat(i).centroid = centroid(i,:);
    mbrStat(i).orientation = dtheta(i);
    idx = (i-1)*4 + 1;
    mbrStat(i).corners = corners(idx:idx+3,:);
end

%MBRvisualization(mbrStat,mov,1)
%%
figure('Position',[164 280 1397 474])
subplot(1,2,1)
h1 = imshow(mov(1).gray,[1 255]);
hold on
h2 = line([centroid(1,1) centroid(1,1)+dxdy(1,1)],[centroid(1,2) centroid(1,2)+dxdy(1,2)],'Color','g');
title('Velocity Vector Frame: 1')
t1 = text(40,900,num2str(bodyTheta(1)),'FontSize',30,'Color','r');
subplot(1,2,2)
title('BodyFrame Velocity')
h3 = line([0 dxdyBody(1,1)],[0 -dxdyBody(1,2)],'Color','g');% flip y axis b/c of image axis
axis([-1400 1400 -1400 1400])
axis equal


if false
for l = 2:169
    subplot(1,2,1)
    title(strcat('Velocity Vector Frame: ',num2str(l)))
    set(h1,'CData',mov(l).gray);
    set(h2,'XData',[centroid(l,1) centroid(l,1)+dxdy(l,1)],'YData',[centroid(l,2) centroid(l,2)+dxdy(l,2)])
    set(t1,'String',num2str(bodyTheta(l)))
    set(h3,'XData',[0 dxdyBody(l,1)],'YData',[0 -dxdyBody(l,2)])
    
    %keyboard
end
end
figure
title('Histogram of Body Frame Velocity Vectors')
load('frame1Rotated.mat');
% shift to 0,0
% centroid = [621,590];
imgOffsetXY = [621,590];

imshow(imgRot,[1 255])
hold on
idxlines = [100:169];
dxscale = 5;
linex = [zeros(1,length(idxlines));dxdyBody(idxlines,1)'/dxscale];
linex = bsxfun(@plus,linex,imgOffsetXY(1));

liney = [zeros(1,length(idxlines));dxdyBody(idxlines,2)'/dxscale];
liney = bsxfun(@plus,liney,imgOffsetXY(2));

line(linex,liney)

%% Generate image of frames 100 - 170
figure('Position',[164 280 1397 474])
subplot(1,2,1)
imshow(mov(100).gray,[1 255])
subplot(1,2,2)
imshow(mov(170).gray,[1 255])

imgCombined = mov(100).gray/2 + mov(170).gray/2;

figure
imshow(imgCombined,[1 255])
hold on
for i = 100:length(centroid)
    plot(centroid(i,1),centroid(i,2),'Color',[(i - 100)/70 (170-i)/70 0],'Marker','.','MarkerSize',15)
    %plot(centroid(100:end,1),centroid(100:end,2),'.b','MarkerSize',15)
end
title('2H_40X.avi Trajectory Frames 100 - 170')

%% print -dtiff TrajectoryFrames100_170
thetaSave = theta(100:170) - theta(100);
H2angle = thetaSave;
H2xy = bsxfun(@minus,centroidMicron(100:170 ,:),centroidMicron(100,:));
save('2Htrackstat','H2angle','H2xy')


