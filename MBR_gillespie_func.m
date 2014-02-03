function [timeVec,MBRstate]= MBR_gillespie_func(rxnrate,init,simIterations,attractant,MBRcorners,varargin)
% Used in EcoliSimMain to run n-independent Gillespie simulations for cells
% adhered to an MBR.
% ----- Inputs
%       - reaction rates, rxnrate - vector a reaction rates in cell
%       - intial state of chemicals in each cell (same)
%       - number os simulation Iterations, simIterations (time will be
%       equalized in post processing when comparing multiple simulations
%       - attractant function that results in change in reaction rates

% ----- Outputs
%       - timeVec representing when the reactions occur
%       - MBRstate vector with information on state of MBR corresponding to the
%       timeVec.

%MBRstate = struct('posn',[0 0 360*rand(1)],'cellposn',[],'F',[]); %in fixed/world frame
%MBRstate = struct('posn',[0 0 90-9],'cellposn',[],'F',[]); %in fixed/world frame
MBRstate = struct('posn',[0 0 60],'cellposn',[],'F',[]); %in fixed/world frame

% varargin for cell distribution provided
if nargin == 6 && ~isempty(varargin{1})
    MBRstate.cellposn = varargin{1};
    
else
    MBRstate.cellposn = [];
    numcell = 200;
end

% tumble torque
%tumbleconstant = 1e-9; %[Nm] CCW>0, CW<0

tumbleconstant = 0; %[Nm] CCW>0, CW<0
%Tumble constant should be positive, so that it exerts CCW torque on MBR
%update MBR state
%p = 0.41e-12; %pN from lit
p = 0.45e-12; %pN from lit
% convert p and q to up and uq
%p = p*1e6;
%q = q*1e6;

% EBS's tranlating U COMSOL
%FdragTrans = 7e-12; %N = kg m/s^2 (7pN)
FdragTrans = 1.50e-11; %N
VdragTrans = 10e-6; %m/s,

kt = FdragTrans/VdragTrans % units Ns/m

%% 
%E EB's rotation COMSOL for sq.
%FdragRot = 2.46e-14; %Nm, torque on bottom, F*d EB comsol 10^-14Nm, 10^-2 pNm
%FdragRot = 2.4e-17; %Nm, DW back of envelope based on translation

%% new model
%FdragRot = 6.03e-13; % N from NEW COMSOL 
FdragRot = 3.29e-13; % N from NEW COMSOL 
%VdragRot = 10; %deg/s
VdragRot = rad2deg(0.016); %deg/s 0.016 rad / s
%kr = FdragRot*30e-6/deg2rad(VdragRot) % kg m^2/s 
kr = FdragRot*30/VdragRot % N um s/deg

%% compute q force based on torque and number of edge cells
q = FdragRot;  % trans

% LEAST SQS P.Q
p = .4804e-12;
q = -.8258e-12;

p =   0.5070e-14;
q =  -0.6292e-14;

p = 0.5701e-8;
q = 0.0941e-8;

%p = .4804e-12;
%q = -.8258e-12;

%%
%uFdragRot = FdragRot * (1e12);
%VdragRot = rad2deg(0.016); %rad/s
%kr = FdragRot/VdragRot; % kg m^2/s 

%Initiate cells on MBR (default)
celllength = 10; %um

%initialize variables
timeVec(1) = 0;
rxncell(1) = 0;

if isempty(MBRstate.cellposn)
    [MBRstate.cellposn,edgecell,bacHead,bacTail] = mbr_cell_distribution(MBRcorners,numcell,celllength);
else
    numcell = length(MBRstate.cellposn);
end
% initial chem state:
cellstate = repmat(init,[1,1,numcell]); %cellstate(timeidx,chem,cell)
% matrix storing the time of next rxn in row 1 and rxn number in row 2
nextRxnTime = zeros(2,numcell);

% randomize start state (run or tumble)
runTumble = rand(1,numcell);
MBRstate.F = ones(1,numcell);
MBRstate.F(runTumble<1/11) = 0;
%keyboard

% convert from um to m:
%MBRcorners.cells = MBRcorners.cells*10^-6;
%MBRcorners.nocells = MBRcorners.nocells*10^-6;
%MBRstate.cellposn(:,1:2) = MBRstate.cellposn(:,1:2)*10^-6;
%celllength = celllength*10^-6;

%% Determine if the bacterium hangs over the edge of microstructure
% determine edge bacteria
%edgecell = zeros(1,numcell); % store 1 if edge bacterium

[edgecell,bacHead,bacTail] = find_edge_bacteria(MBRcorners.cells,MBRcorners.nocells,MBRstate.cellposn,celllength);

drawlength = 3;
[~,~,bacTailDraw] = find_edge_bacteria(MBRcorners.cells,MBRcorners.nocells,MBRstate.cellposn,drawlength);

%% test correct edge detection
plotForceVec = false;

if false
    flagella1 = figure;
    
    % draq sq
    x1 = MBRcorners.cells(1,1);
    y1 = MBRcorners.cells(1,2);
    x2 = MBRcorners.cells(2,1);
    y2 = MBRcorners.cells(2,2);
   
    %cornerallpts(:,1) = [x1,x2,x2,x1,x1];
    %cornerallpts(:,2) = [y1,y1,y2,y2,y1];
    
    % draw h
    xmp = 18; %midpoint for H
    ymp = 8;
    cornerallpts(:,1) = [x1,-xmp, -xmp, xmp, xmp, x2,x2, xmp xmp -xmp -xmp x1,x1];
    cornerallpts(:,2) = [y1, y1, -ymp, -ymp, y1,  y1,y2, y2, ymp, ymp, y2, y2,y1];
    
    %% draw h-Htranslating
    %xmp = 18; %midpoint for H
    %ymp = 9;
    %x2 = 25;
    
    %cornerallpts(:,1) = [x1,-xmp, -xmp, xmp-5, xmp-5, x2,x2, xmp-5 xmp-5 -xmp -xmp x1,x1];
    %cornerallpts(:,2) = [y1, y1, -ymp, -ymp, y1,  y1,y2, y2, ymp, ymp, y2, y2,y1];
    %%
    
    plot(cornerallpts(:,1),cornerallpts(:,2),'-k','LineWidth',4)
    
    for cell = 1:numcell
        % determine if it is in the MBR
        hold on
        bacX = [bacHead(cell,1);bacTailDraw(cell,1)];
        bacY = [bacHead(cell,2);bacTailDraw(cell,2)];
        hold on
        if edgecell(cell)
            % OVER EDGE IN RED
            cellplot(cell) = plot(bacX,bacY,'-r','LineWidth',5);
        else
            cellplot(cell) = plot(bacX,bacY,'-b','LineWidth',5);
        end
    end
    axis equal
    plot(bacHead(:,1),bacHead(:,2),'.k')        
    plot(bacHead(:,1),bacHead(:,2),'ok','MarkerSize',10)    
    if plotForceVec
        % propulsion force
        th = MBRstate.cellposn(:,3);
        dxNormal = -cosd(th);
        dyNormal = -sind(th);
        quiver(MBRstate.cellposn(1:numcell,1), MBRstate.cellposn(1:numcell,2),dxNormal,dyNormal);
        
        % side force
        thTangent = edgecell.* MBRstate.cellposn(:,3);
        dxTangent = sind(thTangent);
        dyTangent = -cosd(thTangent);
        quiver(MBRstate.cellposn(edgecell==1,1), MBRstate.cellposn(edgecell==1,2),dxTangent(edgecell==1),dyTangent(edgecell==1));
    end
    axis ij
    print -dtiff 'CellDist'
    keyboard
    
end

%% generate a time and rxn for each cell
for j = 1:numcell
    %determine what reaction occurs and return new state and change in cell
    %angle and tau time
    [tau,newcellchem,~] = cell_gillespie(timeVec,rxnrate,init,'MBR');
    nextRxnTime(:,j) = [tau;newcellchem];
end

%% trouble shooting dxdt_f, dydt_f, dadt_f, dxdt, dydt, dadt
troubleshoot = false;
if troubleshoot
    disp('Troubleshoot = True')
    figure;
    dxdt_fvec = 0;
    dydt_fvec = 0;
    dadt_fvec = 0;
    
    dxdtvec = 0;
    dydtvec = 0;
    dadtvec = 0;
    
    subplot(2,1,1)
    hforcexf = plot([1],[0],'xb');
    hold on
    hforcexw = plot(1,0,'ob');
    hforceyf = plot(1,0,'xr');
    hforceyf = plot(1,0,'or');
    xlabel('FrameNumber')
    ylabel('Force')
    subplot(2,1,2)
    hforcea = plot(1,0,'xk');
    xlabel('FrameNumber')
    ylabel('Angular Velocity')
    
    axis ij
end

%disp('MBR simulation started')

% initialize cellAngles
MBRstate.cellAngle(1,:) = MBRstate.cellposn(:,3)';


%% cycle through reactions to determine dynamics
for i = 2:simIterations
    % determine next reaction from nextRxnTime
    [rxntime,cellnum] = min(nextRxnTime(1,:));
    rxncell(i) = cellnum;
    timeVec(i) = rxntime;
    
    %update all cell states
    cellstate(i,:,:) = cellstate(i-1,:,:);
    cellstate(i,:,rxncell(i)) = MBRupdatestate(cellstate(i,:,rxncell(i)),nextRxnTime(2,rxncell(i)));
    
    % update nextRxnTime for rxncell(i)
    [tau,newcellchem,~] = cell_gillespie(timeVec,rxnrate,cellstate(i,:,rxncell(i)),'MBR');
    
    nextRxnTime(:,rxncell(i)) = [timeVec(i) + tau;newcellchem];
    
    bx = MBRstate.cellposn(:,1);
    by = MBRstate.cellposn(:,2);
    th = MBRstate.cellposn(:,3);
    %thunwrap = rad2deg(unwrap(deg2rad(th)));
    %th = thunwrap;
    
    Fnow = MBRstate.F(i-1,:)';
    
    %% Force Calc -> dynamics
    % no side force
    %dxdt_f = kt*sum(Fnow.*cosd(th)) ; %mbr frame
    %dydt_f = kt*sum(Fnow.*sind(th)); %mbr frame
    %numtumblecells = sum(Fnow==0);
    %dadt_f = kr*sum(bx.*Fnow.*sind(th) + by.*Fnow.*cosd(th)) - tumbleconstant*numtumblecells;
    
    % with side force
    dxdt_f = 1/kt*(-p*sum(Fnow.*cosd(th)) + q*sum(Fnow.*edgecell.*sind(th))); %mbr frame
    dydt_f = 1/kt*(-p*sum(Fnow.*sind(th)) - q*sum(Fnow.*edgecell.*cosd(th))); %mbr frame
      
    numtumblecells = sum(Fnow==0);
    dadt_f = 1/kr*(sum(bx.*Fnow.*sind(th)*p + by.*Fnow.*cosd(th)*p + q*Fnow.*edgecell.*(-bx.*cosd(th) + by.*sind(th))) + tumbleconstant*numtumblecells);
    
    %fprintf('Before: xbody: %8.8f ybody: %8.8f thbody: %8.6f \n',dxdt_f,dydt_f,dadt_f)
    %keyboard
      
    % in world-fixed frame
    dxdt = (dxdt_f)*cosd(MBRstate.posn(i-1,3)) - dydt_f*sind(MBRstate.posn(i-1,3));
    dydt = (dxdt_f)*sind(MBRstate.posn(i-1,3)) + dydt_f*cosd(MBRstate.posn(i-1,3));
    dadt = dadt_f;
    
    if troubleshoot
        % update velocity vectors
        dxdt_fvec = [dxdt_fvec dxdt_f];
        dydt_fvec = [dydt_fvec dydt_f];
        dadt_fvec = [dadt_fvec dadt_f];
        
        dxdtvec = [dxdtvec dxdt];
        dydtvec = [dydtvec dxdt];
        
        frame = 1:i;
        % update figure handles
        set(hforcexf,'XData',frame,'YData',dxdt_fvec);
        set(hforceyf,'XData',frame,'YData',dydt_fvec);
        set(hforcexf,'XData',frame,'YData',dxdt_fvec);
        set(hforceyf,'XData',frame,'YData',dydt_fvec);
        set(hforcea,'XData',frame,'YData',dadt_fvec);
        keyboard
    end
    
    % update position of MBR
    dt = timeVec(i) - timeVec(i-1);
    MBRstate.posn(i,:) = MBRstate.posn(i-1,:) + dt*[dxdt dydt dadt];
    MBRstate.posn(i,3) = mod(MBRstate.posn(i,3),360);
    
    % update tumbled bacterium
    %MBRstate.cellposn(:,3) = MBRstate.cellposn(:,3);
    %MBRstate.cellposn(rxncell(i),3) = mod(MBRstate.cellposn(rxncell(i),3) + dth,360);
    
    % determine if ligand sensed more stuff
    % location of each cell at i:
    
    % if the cell tumbled at the end of this turn, set force to 0
    MBRstate.F(i,:) = cellstate(i,8,:);
    
    MBRstate.cellAngle(i,:) = MBRstate.cellposn(:,3)';
    
end

% convert m to um
MBRstate.posn(:,1:2) = MBRstate.posn(:,1:2)*10^6;


%% Compute deterministic model and save in MBRstate.detPosn
[MBRx,MBRy,MBRth,timeaxis] = runDeterministicModel(kt,kr,p,q,MBRstate.cellposn,edgecell,MBRstate.posn(1,:));
MBRstate.detTime = timeaxis;
MBRstate.detPosn = [MBRx' MBRy' MBRth'];
