% Solve least squares for translating case

load('HtransdxdyBody.mat')
% dxdyBody
% [ dx/dt dy/dt] in um/s
% bodyTheta (in deg) 

fps = 5;

dtheta = -diff(bodyTheta)*fps; %deg/s flip into 
dxdyBody = dxdyBody * fps; % um/s (not um/frame)

nFrames = length(dtheta);

temp = [dxdyBody dtheta'];
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
kt = 1/1.3e-12; % 1/kt

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

Apq = repmat(Arowspq,[nFrames 1]);

disp('Translation Analysis PQ')
xpq = Apq\B

%%%%%%%%%
Atrans = A;
Btrans = B;
Apqtrans = Apq;

%
disp('runDeterministicModel')

% determine time variables
time = 65000;
timestep = 1/1000;

pbar = 10/11*xpq(1);
qbar = xpq(2);

%Initialize posn
r1 = [0;0;0]; %r(1) = r_x, r(2)=r_y, r(3) = phi
r2 = [0;0;0];
rhistory = zeros(3,time);

%$================== FINISH THIS!!!!

for i = 1:time
    Adeterministic = [Arowspq zeros(3,1)];
    dx = Adeterministic * r1;
    xdot = dx(1);
    ydot = dx(2);
    phidot = dx(3);
    
    xdot = (pbar *B1+qbar*G1)*cos(r1(3)) - (pbar*B2+qbar*G2)*sin(r1(3));
    ydot = (pbar *B1+qbar*G1)*sin(r1(3)) + (pbar*B2+qbar*G2)*cos(r1(3));
    phidot = pbar*B3 + qbar*G3;
    
    r2(1) = r1(1) + (timestep)*xdot;
    r2(2) = r1(2) + (timestep)*ydot;
    r2(3) = r1(3) + (timestep)*phidot;
    
    rhistory(:,i) = r2;
    r1 = r2;
end

timeaxis = timestep*[1:time];

figure
subplot(2,1,1)
plot(rhistory(1,:),rhistory(2,:))
ylabel('y','fontsize',20)
xlabel('x','fontsize',20)


subplot(2,1,2)
plot(timeaxis,rhistory(3,:))
ylabel('\phi','fontsize',20)
xlabel('Time [s]','fontsize',20)


%%
disp('Dont run rotation')
break
%%%%%% ROTATION %%%%%%%%
% run rotation least sqs
HrotleastsqsParamFit
Aboth = [Atrans;A];
Bboth = [Btrans;B];
Apqboth = [Apqtrans;Apq];

disp('Both Analysis x(4)')
x = Aboth\Bboth
error = abs(Bboth - Aboth*x);
totError = sum(error);
disp('Total Mean error')
errorX = mean(error(1:3:end));
errorY = mean(error(2:3:end));
errorDeg = mean(error(3:3:end));

%BOTH
%x =
%
%   -0.0006
%    0.0012
%   -0.0658
%    0.2155

disp('P Q fit both rotation and translation')

xpqBoth = Apqboth\Bboth
% xpq =
% 
%    1.0e-13 *
% 
%     0.2503
%     0.1316

deterministicmodel