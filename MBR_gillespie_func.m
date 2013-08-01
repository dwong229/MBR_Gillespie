function [timeVec,MBRstate]= MBR_gillespie_func(rxnrate,init,simIterations,attractant)
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

% tumble torque
tumbleconstant = 0; % CCW>0, CW<0 
%Tumble constant should be positive, so that it exerts CCW torque on MBR

%Initiate cells on MBR
MBRstate = struct('posn',[0 0 360*rand(1)],'cellposn',[],'F',ones(1,4)); %in fixed/world frame
numcell = 4;

%initialize variables
timeVec(1) = 0;
rxncell(1) = 0;
% matrix storing the time of next rxn in row 1 and rxn number in row 2
nextRxnTime = zeros(2,numcell);


% MBR geometry:sq of sizes 40:
%place cells at center of each edge
MBRstate.cellposn(1:numcell,1:2) = [20,0;0,-20;-20,0;0,20]; %in MBR frame
MBRstate.cellposn(1:numcell,3) = [360*rand(numcell,1)];
%MBRstate.cellposn(1:4,3) = [0 -90 -180 -270];
MBRstate.cellposn(1:numcell,3) = [zeros(1,numcell)];
    
% generate a time and rxn for each cell
for j = 1:numcell
    %determine what reaction occurs and return new state and change in cell
    %angle and tau time
    [tau,newcellchem,~] = cell_gillespie(timeVec,rxnrate,init,'MBR');
    nextRxnTime(:,j) = [tau;newcellchem];    
end

% initial chem state:
cellstate = repmat(init,[1,1,numcell]); %cellstate(timeidx,chem,cell)


% trouble shooting dxdt_f, dydt_f, dadt_f, dxdt, dydt, dadt
troubleshoot = true;
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
    
end
    
MBRstate.cellAngle(1,:) = MBRstate.cellposn(:,3)';

for i = 2:simIterations
    % determine next reaction from nextRxnTime
    
    [rxntime cellnum] = min(nextRxnTime(1,:));
    rxncell(i) = cellnum;
    timeVec(i) = rxntime;
    
    %update all cell states
    cellstate(i,:,:) = cellstate(i-1,:,:);
    cellstate(i,:,rxncell(i)) = MBRupdatestate(cellstate(i,:,rxncell(i)),nextRxnTime(2,rxncell(i)));
    
    % update nextRxnTime for rxncell(i)
    [tau,newcellchem,~] = cell_gillespie(timeVec,rxnrate,cellstate(i,:,rxncell(i)),'MBR');
    
    nextRxnTime(:,rxncell(i)) = [timeVec(i) + tau;newcellchem];    
    
    %update MBR state
    % viscosity constants
    kt = 2; % p/kt
    kr = 1;
    
    bx = MBRstate.cellposn(:,1);
    by = MBRstate.cellposn(:,2);
    th = MBRstate.cellposn(:,3);
    %thunwrap = rad2deg(unwrap(deg2rad(th)));
    %th = thunwrap;
    
    Fnow = MBRstate.F(i-1,:)';
    
    dxdt_f = kt*sum(Fnow.*cosd(th)); %mbr frame
    dydt_f = kt*sum(Fnow.*sind(th)); %mbr frame
    numtumblecells = sum(Fnow==0);

    dadt_f = kr*sum(bx.*Fnow.*sind(th) + by.*Fnow.*cosd(th)) - tumbleconstant*numtumblecells;
    
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
