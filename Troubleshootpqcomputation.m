%% Test model generated by deterministic model to compute p and q using least squares
clear 
load('TestTrajFromDeterministicModel.mat')
% world positon and orientation

% Compute dxdybody
dxdyBody = repmat([-1.0816e-5, -2.0264e-5],length(timeaxis),1); %[-1.0816e-5m/s, -2.0264e-5m/s, -8.3981 deg/s]

fps = 5;
timeVecTrack = timeaxis;
rawdxdyBody = dxdyBody;

scaleVelocity = 1;


dtheta = ones(1,length(timeaxis))*-8.3981; %deg/s flip into
%dxdyBody = dxdyBody * fps * 10e-6* scaleVelocity; % um/s -> m/s (not um/frame)

dxdyBody(:,1) = smooth(dxdyBody(:,1),3,'lowess');
dxdyBody(:,2) = smooth(dxdyBody(:,2),3,'lowess');

nFrames = length(dtheta);


% What we are trying to fit
vh = figure;
subplot(3,1,1)
plot(timeVecTrack,rawdxdyBody(:,1),'.r')
ylabel('X velocity (um/s)')
title('Velocity in the body frame coord [rot]')
subplot(3,1,2)
plot(timeVecTrack,rawdxdyBody(:,2),'.r')
ylabel('Y velocity (um/s)')
subplot(3,1,3)
%plot([timeVecTrack,timeVecTrack(end)+1/fps] ,bodyTheta,'.r')
plot(timeVecTrack,dtheta,'.r')
ylabel('Angular Velocity (deg/s)')
xlabel('time (s)')

keyboard

temp = [dxdyBody*1 dtheta'*1];
tempT = temp';

% Ax = B
B = tempT(:);

% load cellposn
%load('headangle_data_2H_40X.mat')
load('cellposnOpenCV2H_headangle.mat')

%translating H
MBRcorners.cells(:,1) = [-30;25]; %x coordinates
MBRcorners.cells(:,2) = [-30;25]; %y coordinates

MBRcorners.nocells = [-18 8;13 30;...
    -12 -30;13 -5];

celllength = 10;
[edgecell,~,~] = find_edge_bacteria(MBRcorners.cells,MBRcorners.nocells,cellposn,celllength);

th = cellposn(:,3);
bx = cellposn(:,1);
by = cellposn(:,2);

% define A matrix assuming all cells are exerting a force
Arows = zeros(3,4);

% sum cos(th)
Arows(1,1) = -sum(cosd(th)); % negative because propulsion is from tail to head

% sum sin(th)
Arows(1,2) = sum(edgecell.*sind(th));

% sum cos(th)
Arows(2,1) = -sum(sind(th)); % negative because propulsion is from tail to head

% sum sin(th)
Arows(2,2) = -sum(edgecell.*cosd(th));

% sum b*cos(th)
Arows(3,3) = sum(bx.*sind(th) + by.*cosd(th));

% sum b*sin(th)
Arows(3,4) = -sum(edgecell.*(bx.*cosd(th) + by.*sind(th)));

A = repmat(Arows,[nFrames 1]);

x = A\B

% x =
%
%    -0.0006
%     0.0015
%          0
%     0.1863
%

disp('Translation Fit Error analysis')
error = abs(B - A*x);
% totError = sum(error)
% errorX = sum(error(1:3:end))
% errorY = sum(error(2:3:end))
% errorDeg = sum(error(3:3:end))

%%%%% solve for p.q.

kr = 1/1.0766e-11; % 1/kr
%kt = 1/1.3e-12; % 1/kt
kt = 1/1.5e-6; % 1/kt

% compute parameters
B1 = - kt * sum(cosd(th));
B2 = -kt * sum(sind(th));
B3 = kr * sum(bx.*sind(th) + by.*cosd(th));
G1 = kt * sum(edgecell.*sind(th));
G2 = - kt * sum(edgecell.*cosd(th));
G3 = kr * sum(-edgecell.*bx.*cosd(th) + edgecell.*by.*sind(th));

Arowspq = zeros(3,2);

% sum cos(th)
Arowspq(1,1) = -kt*sum(cosd(th)); % negative because propulsion is from tail to head

% sum sin(th)
Arowspq(1,2) = kt*sum(edgecell.*sind(th));

% sum cos(th)
Arowspq(2,1) = -kt*sum(sind(th)); % negative because propulsion is from tail to head

% sum sin(th)
Arowspq(2,2) = -kt*sum(edgecell.*cosd(th));

% sum b*cos(th)
Arowspq(3,1) = kr*sum(bx.*sind(th) + by.*cosd(th));

% sum b*sin(th)
Arowspq(3,2) = kr*sum(edgecell.*(-bx.*cosd(th) + by.*sind(th)));

% scale velocity:
Arowspq(1:2,:) = Arowspq(1:2,:)*scaleVelocity;

Apq = repmat(Arowspq,[nFrames 1]);


disp('Translation Analysis PQ')
xpqCompute = Apq\B;
% multiply by 11/10 to correct for run:tumble
xpqCompute = xpqCompute*11/10;


