function [MBRx,MBRy,MBRth,timeaxis] = runDeterministicModel(kt_input,kr_input,runforce,sideforce,cellposn,edgecell,initialposn)

% inputs

% outputs

% Take reciprocal.
kt = 1/kt_input;
kr = 1/kr_input;

plotOption = false;

%Simulation of Serratia Marcescens blotted on SU8 rectangular substrate
%Based on Modeling, control and experimental characterization of
%microbiorobots - Selman Sakar et al. (2011)

%System : 50um x 100um, B1 = 13.03um/(s pN), B2 = -43.64um/(s pN),
%B3 = 1.24 (rad/s pN)
% determine time variables
time = 10000;
timestep = 1/1000;

% unpack inputs
pbar = 10/11*runforce;
qbar = 10/11*sideforce;

th = cellposn(:,3);
bx = cellposn(:,1);
by = cellposn(:,2);

%Initialize posn
r1 = [initialposn']; %r(1) = r_x, r(2)=r_y, r(3) = phi
r2 = [initialposn'];
rhistory = zeros(3,time);

% compute parameters
B1 = - kt * sum(cosd(th));
B2 = -kt * sum(sind(th));
B3 = kr * sum(bx.*sind(th) + by.*cosd(th));
G1 = kt * sum(edgecell.*sind(th));
G2 = - kt * sum(edgecell.*cosd(th));
G3 = kr * sum(-edgecell.*bx.*cosd(th) + edgecell.*by.*sind(th));


%% Populate A matrix
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
Arowspq(3,2) = -kr*sum(edgecell.*(bx.*cosd(th) + by.*sind(th)));

Adeterministic = [Arowspq zeros(3,1)];

for i = 1:time
    
    %dx = Adeterministic * r1;
    %xdot = dx(1);
    %ydot = dx(2);
    %phidot = dx(3);
       
    %% 4:46pm 2/2/2014 changes
    %G1 = kt * sum(sind(th));
    xbody = (pbar *B1+qbar*G1);
    ybody = (pbar*B2+qbar*G2);
    
    xdot = (pbar *B1+qbar*G1)*cosd(r1(3)) - (pbar*B2+qbar*G2)*sind(r1(3));
    ydot = (pbar *B1+qbar*G1)*sind(r1(3)) + (pbar*B2+qbar*G2)*cosd(r1(3));
    phidot = pbar*B3 + qbar*G3;
    %disp('Deterministic')
    %fprintf('xbody: %8.8f ybody: %8.8f phibody: %8.8f \n',xbody,ybody,phidot)
    %keyboard
    
    r2(1) = r1(1) + (timestep)*xdot;
    r2(2) = r1(2) + (timestep)*ydot;
    r2(3) = r1(3) + (timestep)*phidot;
    
    rhistory(:,i) = r2;
    r1 = r2;
end

timeaxis = timestep*[1:time];


if plotOption
       
    figure
    subplot(2,1,1)
    plot(rhistory(1,:),rhistory(2,:))
    ylabel('y','fontsize',20)
    xlabel('x','fontsize',20)
    hold on
    plot(rhistory(1,1),rhistory(2,1),'xg')
    
    
    subplot(2,1,2)
    plot(timeaxis,rhistory(3,:))
    ylabel('\phi','fontsize',20)
    xlabel('Time [s]','fontsize',20)
    
end

MBRx = rhistory(1,:);
MBRy = rhistory(2,:);
MBRth = rhistory(3,:);
