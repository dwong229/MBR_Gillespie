% Solve least squares for H rotating case

load('HRotdxdyBody.mat')
% dxdyBody
% [ dx/dt dy/dt] in um/s
% bodyTheta (in deg) 

fps = 4;

dxdyBody = dxdyBody(30:79,:);  %take last 50frames

dtheta = diff(bodyTheta)*fps; %deg/s
dtheta = -dtheta(30:79);
dxdyBody = dxdyBody*fps;

kr = 1/1.0766e-11; % 1/kr
kt = 1/1.5e-6; % 1/kt

nFrames = length(dtheta);

temp = [dxdyBody dtheta'];
tempT = temp';

% Ax = B
B = tempT(:);

% load cellposn
load('headangle_data_H3reverse.mat')

  MBRcorners.cells(:,1) = [-30;30]; %x coordinates
        MBRcorners.cells(:,2) = [-30;30]; %y coordinates
        
        MBRcorners.nocells = [-18 8;18 30;-18 -30;18 -8];        
        
        
celllength = 10;
[edgecell,~,~] = find_edge_bacteria(MBRcorners.cells,MBRcorners.nocells,cellposn,celllength);


th = cellposn(:,3);
bx = cellposn(:,1);
by = cellposn(:,2);
% define A matrix assuming all cells are exerting a force
Arows = zeros(3,4);

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

A = repmat(Arowspq,[nFrames 1]);

disp('Rotation Analysis x(4)')
x = A\B
% x =
% 
%     0.0071
%    -0.0349
%    -0.0035
%          0

% compute error
disp('Error analysis')
error = abs(B - A*x);
totError = sum(error)
errorX = sum(error(1:3:end))
errorY = sum(error(2:3:end))
errorDeg = sum(error(3:3:end))

%%%%% solve for p.q.

