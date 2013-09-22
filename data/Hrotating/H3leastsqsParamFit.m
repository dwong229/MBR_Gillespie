% Solve least squares for H3 translating case

load('H3dxdyBody.mat')
% dxdyBody
% [ dx/dt dy/dt] in um/s
% bodyTheta (in deg) 

fps = 10;

dtheta = diff(bodyTheta)*fps; %deg/s

temp = [dxdyBody dtheta'];
tempT = temp';

% Ax = B
B = tempT(:);

% load cellposn
load('.mat')

% define A matrix assuming all cells are exerting a force
Arows = zeros(3,4);

% sum cos(th)

% sum sin(th)

% sum b*cos(th)

% sum b*sin(th)
