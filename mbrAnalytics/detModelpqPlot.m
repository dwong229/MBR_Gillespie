% Run simulation for multiple p,q and dtermine
close all

global invkt invkr cellposn edgecell pmin pmax qmin qmax p q
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



% rotating H

celllength = 10;
[edgecell,~,~] = find_edge_bacteria(MBRcorners.cells,MBRcorners.nocells,cellposn,celllength);


%% Input force range to test
pmin = 1e-15;
pmax = 1e-12;
pres = 10;

qmin = -1e-13;
qmax = 1e-13;
qres = 10;

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

plotIdx = 0;
for pIdx = 1:length(pVec)
    p = pVec(pIdx);
    for qIdx = 1:length(qVec)
        plotIdx = plotIdx + 1;
        q = qVec(qIdx);
        
        % run deterministic sim
        [~,~,~,~,dxdt_body] = runDeterministicModel(1/invkt,1/invkr,p,q,cellposn,edgecell,[0;0;0]');
        % returns in m
        
        % unpack velocities
        dxdt(plotIdx) = dxdt_body(1) * 1e6;
        dydt(plotIdx) = dxdt_body(2) * 1e6;
        dthdt(plotIdx) = dxdt_body(3);
        
        % assign plot axes
        pplot(plotIdx) = p;
        qplot(plotIdx) = q;
    end
end

StartIdx = floor(length(pplot)/2);

h1 = figure('Position',[113 302 1484 505]);
subplot(1,2,1)
plot3(pplot,qplot,dxdt,'.r')
hold on
plot3(pplot,qplot,dydt,'.b')
dxcircle = plot3(pplot(StartIdx),qplot(StartIdx),dxdt(StartIdx),'or');
dycircle = plot3(pplot(StartIdx),qplot(StartIdx),dydt(StartIdx),'ob');
legend('dxdt','dydt','Location','NorthWest')
axis square
title(strcat(plotstr,': linear velocity'))
xlabel('p-force (N)')
ylabel('q-force (N)')
zlabel('velocity (um/s)')

subplot(1,2,2)
plot3(pplot,qplot,dthdt,'.b')
% Add plot gui
hold on
dthcircle = plot3(pplot(StartIdx),qplot(StartIdx),dthdt(StartIdx),'oc');
axis square
title('angular velocity')
xlabel('p-force (N)')
ylabel('q-force (N)')
zlabel('angular velcity (deg/s)')
p = pplot(floor(length(pplot)/2));
q = qplot(floor(length(pplot)/2));
%% slider!
uicontrol('Style', 'slider',...
        'Min',0,'Max',100,'Value',50,...
        'Position', [650 420 240 20],...
        'Callback', {@PlotSlider_Callback,dxcircle,dycircle,dthcircle});
    
uicontrol('Style','text',...
        'Position',[650 445 240 20],...
        'String','Select force p')

uicontrol('Style', 'slider',...
        'Min',0,'Max',100,'Value',50,...
        'Position', [650 365 240 20],...
        'Callback', {@QPlotSlider_Callback,dxcircle,dycircle,dthcircle});
uicontrol('Style','text',...
        'Position',[650 390 240 20],...
        'String','Select force q')