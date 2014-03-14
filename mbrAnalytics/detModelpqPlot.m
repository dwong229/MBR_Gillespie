% Run simulation for multiple p,q and dtermine
close all
clear all
global invkt invkr cellposn edgecell pmin pmax qmin qmax p q krmin krmax ktmin ktmax pVec qVec vvector thvector
%% cell dist option
mbrcelloption = 1;

%% Load MBR cell info
% drag coeff
invkr = 1/1.0766e-11; % 1/kr
invkt = 1/1.5e-6; % 1/kt


%translating H
if mbrcelloption == 1;
    disp('Run Translating')
    load('cellposnOpenCV2H_headangle.mat')
    MBRcorners.cells(:,1) = [-30;25]; %x coordinates
    MBRcorners.cells(:,2) = [-30;25]; %y coordinates
    MBRcorners.nocells = [-18 8;13 30;...
        -12 -30;13 -5];
    plotstr = 'H Translating';
    
else
    disp('Run Rotating')
    load('headangle_data_H3reverse.mat')
    MBRcorners.cells(:,1) = [-30;30]; %x coordinates
    MBRcorners.cells(:,2) = [-30;30]; %y coordinates
    MBRcorners.nocells = [-18 8;18 30;-18 -30;18 -8];
    plotstr = 'H Rotating';
end

celllength = 10;
[edgecell,~,~] = find_edge_bacteria(MBRcorners.cells,MBRcorners.nocells,cellposn,celllength);


%% Input force range to test
pmin = 1e-15;
pmax = 1e-12;
pres = 10;

qmin = -1e-13;
qmax = 1e-13;
qres = 10;

krmin = 0; %0.1*(1/invkr);
krmax = 2*(1/invkr);
ktmin = 0; %0.1*(1/invkr);
ktmax = 2*(1/invkt);

%% set up vectors of p and q force to cycle through
pdiff = pmax - pmin;
qdiff = qmax - qmin;

pVec = pmin:pdiff/pres:pmax;
qVec = qmin:qdiff/pres:qmax;

% initialize vectors to store plot info.
pplot = zeros(1,pres*qres);
qplot = zeros(1,pres*qres);
dthdt = zeros(1,pres*qres);
dxdt = zeros(1,pres*qres);
dydt = zeros(1,pres*qres);

% compute velocities
[dxdt,dydt,dthdt,pplot,qplot] = UpdateVelcityDetModel(pVec,qVec);

%% Make initial plots
StartIdx = floor(length(pplot)/2);

h1 = figure('Position',[113 302 1526 505]);
subplot(1,2,1)
xvelplot = plot3(pplot,qplot,dxdt,'.r');
hold on
yvelplot = plot3(pplot,qplot,dydt,'.b');
dxcircle = plot3(pplot(StartIdx),qplot(StartIdx),dxdt(StartIdx),'xk','MarkerSize',8);
dycircle = plot3(pplot(StartIdx),qplot(StartIdx),dydt(StartIdx),'xk','MarkerSize',8);
legend('dxdt','dydt','Location','NorthWest')
axis square
title(strcat(plotstr,': linear velocity'))
xlabel('p-force (N)')
ylabel('q-force (N)')
zlabel('velocity (um/s)')

subplot(1,2,2)
thvelplot = plot3(pplot,qplot,dthdt,'.b');
% Add plot gui
hold on
dthcircle = plot3(pplot(StartIdx),qplot(StartIdx),dthdt(StartIdx),'xk','MarkerSize',8);
axis square
title('angular velocity')
xlabel('p-force (N)')
ylabel('q-force (N)')
zlabel('angular velcity (deg/s)')
p = pplot(floor(length(pplot)/2));
q = qplot(floor(length(pplot)/2));

%% Make MBR plot

mbrplot = figure('Position',[1267,2,560,420])
plot([-20,20,20,-20,-20],[-20,-20,20,20,-20],'-k') % sq robot
axis([-40 40 -40 40])
axis equal
xlabel('X')
ylabel('Y')
title('Velocity Vectors')
vx = get(dxcircle,'ZData');
vy = get(dycircle,'ZData');
th = get(dthcircle,'ZData');

%normv = [vx vy]/norm([vx vy]);
normv = [vx vy]/10;

hold on
%arrowscale = 15;
vvector = quiver(0,0,normv(1),normv(2),'MarkerSize',10,'MaxHeadSize',1);
thvector = quiver(22,0,0,th,'MarkerSize',10,'MaxHeadSize',10);

figure(h1)
%% p and q force text
ptext = uicontrol('Style','text',...
        'Position',[650 445 240 20],...
        'String',strcat('Select force p: ',num2str(p*1e-12),'pN'));
qtext = uicontrol('Style','text',...
        'Position',[650 390 240 20],...
        'String',strcat('Select force q: ',num2str(q*1e-12),'pN'));
    
%% FORCE slider!    
uicontrol('Style', 'slider',...
        'Min',0,'Max',100,'Value',50,...
        'Position', [650 420 240 20],...
        'Callback', {@PlotSlider_Callback,dxcircle,dycircle,dthcircle,ptext});
    
uicontrol('Style', 'slider',...
        'Min',0,'Max',100,'Value',50,...
        'Position', [650 365 240 20],...
        'Callback', {@QPlotSlider_Callback,dxcircle,dycircle,dthcircle,qtext});
    
%% drag

% update kr + kt
%invkr = 1/1.0766e-11; % 1/kr
%invkt = 1/1.5e-6; % 1/kt
%% DRAG sliders! 
%% krkt text
krtext = uicontrol('Style','text',...
        'Position',[650 325 240 20],...
        'String',strcat('Select drag inv(kr): ',num2str(invkt,6),'m/Ns'));
kttext = uicontrol('Style','text',...
        'Position',[650 270 240 20],...
        'String',strcat('Select drag inv(kt): ',num2str(invkr,6),'m/Ns'));

uicontrol('Style', 'slider',...
        'Min',0,'Max',100,'Value',50,...
        'Position', [650 310 240 20],...
        'Callback', {@dragSlider_Callback,'r',krtext,xvelplot,yvelplot,thvelplot,dxcircle,dycircle,dthcircle});
    
uicontrol('Style', 'slider',...
        'Min',0,'Max',100,'Value',50,...
        'Position', [650 255 240 20],...
        'Callback', {@dragSlider_Callback,'t',kttext,xvelplot,yvelplot,thvelplot,dxcircle,dycircle,dthcircle});