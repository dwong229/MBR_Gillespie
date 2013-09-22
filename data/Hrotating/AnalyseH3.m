% analyse H3_40X.avi from file full of corners
clear
close all
fps = 4;

% extract file of corners in pixels
load('CornersH3_frames_120_200.mat')

nFrames = 81;
% frames 120 - 200 analysed

for i = 1:80
    idx = (i-1)*4 + 1;
    
    % find mean by computing centroid of 4 corners
    % [xmean ymean]
    centroid(i,:) = mean(cornersMat(idx:idx+3,:),1);
    
    % find theta by computing rotation matrix between corners of current
    % and last frame
    
    %given inliers, compute r and t
    
    if i > 1
        
        % put current 
        pNow0 = [cornersMat(idx:idx+3,:) ones(4,1)]';
        
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
    pLast0 = [cornersMat(idx:idx+3,:) ones(4,1)]';
    uLast = [centroid(i,:) 1]';
    
    for j = 1:4
        % organize corners so that they correspond
        % take current corner and set to last
        lastcorner(j,:) = cornersMat(idx+j-1,:);
        
        % find next corner that is closest to lastcorner
        distMat = dist([lastcorner(j,:);cornersMat(idx+4:idx+7,:)]');
        [a minIdx] = min(distMat(1,2:end));
        cornersMatTemp(j,:) = cornersMat(idx+3+minIdx,:);
        
    end
    cornersMat(idx+4:idx+7,:) = cornersMatTemp;
    lastcorner = [];
    cornersMatTemp = [];
end

timeVec = [1:i]*1/fps;
subplot(2,1,1)
theta = cumsum(dtheta);
plot(timeVec,theta,'.b')
title('World Frame Angle')
xlabel('Time (sec)')
ylabel('Angle (deg)')

% convert pixels to microns
% 10x Zeiss
% corners: [92 58] [20 121]
%pix2mic = 95.67/60; % pixels per 60um

% mov(5)
% corners: [649 579] [648 682]
pix2mic = 103/60; % pixels per 60um 

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
pts(1,:) = [530 829]; % bottom left H
pts(2,:) = [470 759];% top left H

initialBodyAngle = rad2deg(atan2(abs(diff(pts(:,2))),diff(pts(:,1)))) - 90;

% Add initial body angle to (Angle to bottom of H from horizontal)
bodyTheta = theta + initialBodyAngle;

dxdy= diff(centroidMicron);

for k = 1:length(dxdy)
    R = @(x)[cosd(x) -sind(x);sind(x) cosd(x)]; % Body = R*World
    dxdyBody(k,:) = R(bodyTheta(k))*dxdy(k,:)';
end

dxdyBody = dxdyBody;

% figure
% title('Body Frame Velocity')
% hold on
% plot(timeVec(2:end),dxdyBody(:,1),'-b');
% plot(timeVec(2:end),dxdyBody(:,2),'-g');
% legend('X','Y')
% xlabel('Time (sec)')
% ylabel('Velocity (um/s)')

load('H3.mat')
mbrStat = struct('centroid',zeros(1,2,'double'),'orientation',zeros(1,'double'),'corners',zeros(4,2,'double'));
for i = 1:80
    mbrStat(i).centroid = centroid(i,:);
    mbrStat(i).orientation = dtheta(i);
    idx = (i-1)*4 + 1;
    mbrStat(i).corners = cornersMat(idx:idx+3,:);
end

%MBRvisualization(mbrStat,mov,1)
%%

% figure('Position',[164 280 1397 474])
% subplot(1,2,1)
% h1 = imshow(mov(1).gray,[1 255]);
% hold on
% h2 = line([centroid(1,1) centroid(1,1)+dxdy(1,1)],[centroid(1,2) centroid(1,2)+dxdy(1,2)],'Color','g');
% title('Velocity Vector Frame: 1')
% t1 = text(40,900,num2str(bodyTheta(1)),'FontSize',30,'Color','r');
% subplot(1,2,2)
% title('BodyFrame Velocity')
% h3 = line([0 dxdyBody(1,1)],[0 -dxdyBody(1,2)],'Color','g');% flip y axis b/c of image axis
% axis([-1400 1400 -1400 1400])
% axis equal
% 
% 
% if false
% for l = 120:200
%     subplot(1,2,1)
%     title(strcat('Velocity Vector Frame: ',num2str(l)))
%     set(h1,'CData',mov(l).gray);
%     set(h2,'XData',[centroid(l,1) centroid(l,1)+dxdy(l,1)],'YData',[centroid(l,2) centroid(l,2)+dxdy(l,2)])
%     set(t1,'String',num2str(bodyTheta(l)))
%     set(h3,'XData',[0 dxdyBody(l,1)],'YData',[0 -dxdyBody(l,2)])
%     
%     keyboard
% end
% end

%% Body Frame Velocity Vector
%figure

%% Generate image of frames 
%figure('Position',[164 280 1397 474])
subplot(1,2,1)
startFrame = 150;
title('Histogram of Body Frame Velocity Vectors')
%load('frame1Rotated.mat');
%shift to 0,0
centroid = [621,590];
imgOffsetXY = [621,590];

%imshow(imgRot,[1 255])
hold on
idxlines = [30:79];
dxscale = 5;
linex = [zeros(1,length(idxlines));dxdyBody(idxlines,1)'/dxscale];
linex = bsxfun(@plus,linex,imgOffsetXY(1));

liney = [zeros(1,length(idxlines));dxdyBody(idxlines,2)'/dxscale];
liney = bsxfun(@plus,liney,imgOffsetXY(2));

line(linex,liney)
im1 = mov(startFrame).gray(640:890,450:700);
im2 = mov(200).gray(640:890,450:700);
imshow(im1,[1 255])
subplot(1,2,2)
imshow(im2,[1 255])


imgCombined = im1/2 + im2/2;

figure
imshow(imgCombined,[1 255])
hold on
for i = startFrame-120:length(centroid)
    plot(centroid(i,1),centroid(i,2),'Color',[(i)/81 (81-i)/81 0],'Marker','.','MarkerSize',15)
    %plot(centroid(100:end,1),centroid(100:end,2),'.b','MarkerSize',15)
end
axis ij
axis on
title('H3.avi Trajectory Frames 130 - 200')
print -dtiff TrajectoryFrames130_170xyaxes

%% save datafile with x-y in microns and angle in degrees from 90
thetaSave = theta(startFrame-120:end) - theta(startFrame-120);
H3angle = thetaSave + (+360)*[thetaSave<-360];
H3xy = bsxfun(@minus,centroidMicron(startFrame-120:end ,:),centroidMicron(startFrame-120,:));
save('H3trackstat','H3angle','H3xy')



