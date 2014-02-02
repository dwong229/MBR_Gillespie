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
%kt = 1/1.3e-12; % 1/kt
kt = 1/1.5e-6; % 1/kt

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

[MBRx,MBRy,MBRth,timeaxis] = runDeterministicModel(1/kt,1/kr,xpq(1),xpq(2),cellposn,edgecell);

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


