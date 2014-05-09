% E.coli Simulation Main
clear
close all
%%%%%%%%
simMode = 5; %1: one simulation, 2: repetition
% 6: repetition of 4
%%%%%%%%
repeatSim = 20;

simIterations = 5000;%5000;

%% simulation parameters
delta = 0; % A -> I   reduced ligand detected
alpha = 0.11; % AA -> AAp
beta = 1e8;    % AAp + YY -> AA + YYp autophosphorelation of intracellular kinase (CheA) + phosphorelation of CheY
gamma = 0.031; % YYp -> YY
%gamma = 100; % YYp -> YY
%mu = 0;
mu = 0.015;      % YYp -> motor
rho = 0.015/10;     % run rate (run to tumble)
%rho = (alpha+beta+gamma+mu)/10;

mu = 10; % tumble to run
rho = 1; % run to tumble

% fudge factors that control impact of conc detection:
phi = gamma;   % I-> A   ligand binding
alpha_a = alpha*1.5; % A + AA -> A + AAp
gamma_a = gamma*1.5; % A + YYp -> A + YY
rxnrate = [phi delta alpha beta gamma mu rho alpha_a gamma_a];

% I = 1; % no chemoattractant
% A = 0;
% AA = 1;  %CheA
% AAp = 0; %CheA-p
% YY = 1;  %CheY
% YYp = 0; %CheY-p
% Mot = 1; %CheY-p attached to motor
% Run = 1;

I = 1; % no chemoattractant
A = 0;
%AA = 6700;  %CheA
AA = 0;
AAp = 0; %CheA-p
%YY = 8200;  %CheY
YY = 1;
YYp = 1; %CheY-p
Mot = 1; %CheY-p attached to motor
Run = 1;


x = 0;
y = 0;
theta = 360*rand(1);

init.chem = [I A AA AAp YY YYp Mot Run];
init.dyn = [x y theta]; %start at [x,y,theta]

%simTime = 300;

switch simMode
    
    case 1
        % Computation for single cell:
        disp('Running Simulation Once')
        % Simulate
        
        attractant = @(x) 0;
        [timeVec,state] =  ecoli_gillespie_func(rxnrate,init,simIterations,attractant);
        simTime = floor(timeVec(end)*0.90);
        
        %plot simulation
        plotStates(timeVec,state)
        
    case 2
        disp('Running Repetition of Simulation')
        
        t = CTimeleft(repeatSim);
        for rep = 1:repeatSim
            t.timeleft();
            plotnow = false;
            
            % run simulation
            attractant = @(x) 0;
            [timeVec,state] =  ecoli_gillespie_func(rxnrate,init,simIterations,attractant);
            
            % compute stats
            if rep == 1
                simTime = floor(timeVec(end)*0.75);
            end
            
            [~,~,rtratio,dx] = runtumble(timeVec,state,simTime);
            
            stats.dx(rep,:) = dx;
            stats.rtratio(rep) = rtratio;
        end
        save('randomwalkstats.mat','stats')
        
        % plot statistics
        figure
        hist(stats.rtratio)
        title('Histogram of run/tumble ratio')
        xlabel('run/tumble ratio')
        ylabel('frequency')
        
        figure
        for i = 1:size(stats.dx,1)
            hold on
            line([0 stats.dx(i,1)],[0 stats.dx(i,2)])
        end
        title('Random Walk Planar Displacement')
        xlabel('x (um)')
        ylabel('y (um)')
        
        avgrtratio = mean(stats.rtratio);
        stdrtratio = std(stats.rtratio);
        fprintf('Mean run-tumble ratio: %4.4f, stdev: %4.4f \n',avgrtratio,stdrtratio)
        
        avgdx = mean(stats.dx,1);
        stddx = std(stats.dx,1);
        fprintf('Mean displacement ratio: (%4.0f,%4.0f), stdev: (%4.0f,%4.0f) \n',avgdx(1),avgdx(2),stddx(1),stddx(2))
        
    case 3
        % Computation for single cell:
        disp('Running Simulation Once with Attractant')
        % Simulate
        
        attractant = @(x) attract_exp(x);
        attractant = @(x) attract_lin(x);
        [timeVec,state] =  ecoli_gillespie_func(rxnrate,init,simIterations,attractant);
        simTime = floor(timeVec(end)*0.75);
        
        %plot simulation
        plotStates(timeVec,state)
        
        pathhandle = gcf;
        xmin = min(state.dyn(:,1));
        xmax = max(state.dyn(:,1));
        ymin = min(state.dyn(:,2));
        ymax = max(state.dyn(:,2));
        
        %plot_attract_exp(pathhandle,[xmin xmax ymin ymax])
        plot_attract_lin(pathhandle,[xmin xmax ymin ymax])
        
    case 4
        disp('Running Repetition of Simulation with attractant')
        
        t = CTimeleft(repeatSim);
        for rep = 1:repeatSim
            t.timeleft();
            plotnow = false;
            
            % run simulation
            
            attractant = @(x) attract_exp(x);
            attractant = @(x) attract_lin(x);
            [timeVec,state] =  ecoli_gillespie_func(rxnrate,init,simIterations,attractant);
            
            % compute stats
            if rep == 1
                simTime = floor(timeVec(end)*0.75);
            end
            
            [~,~,rtratio,dx] = runtumble(timeVec,state,simTime);
            
            stats.dx(rep,:) = dx;
            stats.rtratio(rep) = rtratio;
        end
        save('attractexpstats.mat','stats')
        
        % plot statistics
        figure
        hist(stats.rtratio)
        title('Histogram of run/tumble ratio')
        xlabel('run/tumble ratio')
        ylabel('frequency')
        
        figure
        for i = 1:size(stats.dx,1)
            hold on
            line([0 stats.dx(i,1)],[0 stats.dx(i,2)])
        end
        title('Random Walk Planar Displacement')
        xlabel('x (um)')
        ylabel('y (um)')
        
        avgrtratio = mean(stats.rtratio);
        stdrtratio = std(stats.rtratio);
        fprintf('Mean run-tumble ratio: %4.4f, stdev: %4.4f \n',avgrtratio,stdrtratio)
        
        avgdx = mean(stats.dx,1);
        stddx = std(stats.dx,1);
        fprintf('Mean displacement ratio: (%4.0f,%4.0f), stdev: (%4.0f,%4.0f) \n',avgdx(1),avgdx(2),stddx(1),stddx(2))
        
        
    case 5
        disp('Run MBR simulation once')
        % Computation for single cell:
        % Simulate
        
        %disp('Running Simulation Once with Attractant')
        %attractant = @(x) attract_exp(x);
        disp('Running Simulation Once without Attractant')
        attractant = @(x) 0;
        
        init = [I A AA AAp YY YYp Mot Run];
        
        %Initiate cells on MBR
        numcell = 200;
        celllength = 10; %um
        MBRcorners = struct('cells',zeros(2,2),'nocells',[]);
        
        %% U 40x 40
        MBRcorners.cells(:,1) = [-20;20]; %x coordinates
        MBRcorners.cells(:,2) = [-20;20]; %y coordinates
        MBRcorners.nocells = [-10 -10;10 20];
        
        %% Square 40 x 40
        MBRcorners.cells(:,1) = [-20;20]; %x coordinates
        MBRcorners.cells(:,2) = [-20;20]; %y coordinates
        MBRcorners.nocells = [];
        
        %% H 60 x 60
        % rotating H
        MBRcorners.cells(:,1) = [-30;30]; %x coordinates
        MBRcorners.cells(:,2) = [-30;30]; %y coordinates
        
        MBRcorners.nocells = [-18 8;18 30;-18 -30;18 -8];        
        
        %translating H
        MBRcorners.cells(:,1) = [-30;25]; %x coordinates
        MBRcorners.cells(:,2) = [-30;25]; %y coordinates

        MBRcorners.nocells = [-18 8;13 30;...
          -12 -30;13 -5];
        
        %%  Meters to micrometers
        %MBRcorners.cells = MBRcorners.cells*10^-6;
        %MBRcorners.nocells = MBRcorners.nocells*10^-6;
        
        %%
        %cellposnfile = 'cellposn400cells.mat';
        %cellposnfile = 'cellposnU_200cells_trans.mat';
        %cellposnfile = 'cellposn0angle.mat';
        %cellposnfile = 'cellposnborder.mat';
        
        %cellposnfile = 'headangle_data_2H_40X.mat';
        %cellposnfile = 'headangle_data_H3.mat';
        %cellposnfile = 'headangle_data_H3reverse.mat';
        cellposnfile = 'cellposnOpenCV2H_headangle.mat';
        
        %cellposnfile = [];
        if exist(cellposnfile,'file') == 2
            disp('Loading cell-position file')
            load(cellposnfile)
        else
            disp('Generating cell-positions')
            cellposn = mbr_cell_distribution(MBRcorners,numcell,celllength,1);
        end
        
        [timeVec,state] =  MBR_gillespie_func(rxnrate,init,simIterations,attractant,MBRcorners,cellposn);
        simTime = floor(timeVec(end)*0.75);
        
        %plot simulation
        x = state.posn(:,1);
        y = state.posn(:,2);
        th = state.posn(:,3);
        
        hxyth = figure;
        subplot(2,1,1)
        plot(timeVec,x,'-b',timeVec,y,'-r')
        title('MBR simulation')
        xlabel('Time (s)')
        ylabel('Position (um)')
        legend('x','y')
        xlim([0 ceil(timeVec(end))])
        subplot(2,1,2)
        plot(timeVec,th)
        xlabel('Time (s)')
        ylabel('Orientation (deg')
        axis([0,ceil(timeVec(end)),0,360])
        
        % save figure for position and orientation
        print(gcf,'-djpeg','posOrientationSim1.jpg')
        
        htraj = figure;
        plot(x,y,'-g')
        xlabel('x')
        ylabel('y')
        title('Position (um)')
        
        %---- Compare with deterministic model ----
        keyboard
        EndTimeIdx = find(timeVec(end)>state.detTime,1,'last') + 1;
        %timeaxis = state.detTime(1:EndTimeIdx);
        timeaxis = timeVec(1:EndTimeIdx);
        xDet = state.detPosn(1:EndTimeIdx,1)*10^6;
        yDet = state.detPosn(1:EndTimeIdx,2)*10^6;
        thDet = state.detPosn(1:EndTimeIdx,3);
        
        figure(hxyth) 
        subplot(2,1,1)
        hold on
        plot(timeaxis,xDet,'.b')
        plot(timeaxis,yDet,'.r')
        subplot(2,1,2)
        hold on
        plot(timeaxis,thDet,'.b')

        
        figure(htraj)
        hold on
        plot(xDet,yDet,'.b')
        axis ij
        % -- end determinent comparison -- 
        
        
        [rtRatio,dtheta,dx] = eval_MBR_gillespie(timeVec,state,simTime)
        % make a movie to simulate mbr state:
        MBRmovie(timeVec,state)
        
        % run a expt compare file
        %H3compare
        HtransCompare(x,y,th,timeVec,'sim')
        
        
        
    case 6
        disp('Running Repetition of MBR Simulation')
        
        attractant = @(x) 0;
        
        fixedcellposn = true;
        
        if fixedcellposn
            disp('Generate a single distribution and observe variability')
            %Initiate cells on MBR
            numcell = 200;
            celllength = 3; %um
            %MBRcorners.cells = [-20 -20;20 20];
            %MBRcorners.nocells = [];
            
            %% H 60 x 60
            %MBRcorners.cells(:,1) = [-30;30]; %x coordinates
            %MBRcorners.cells(:,2) = [-30;30]; %y coordinates
            
            %MBRcorners.nocells = [-10 10;10 30;-10 -30;10 -10];
            MBRcorners.nocells = [-18 8;18 30;-18 -30;18 -8];        
            
            %translating H
            MBRcorners.cells(:,1) = [-30;30]; %x coordinates
            MBRcorners.cells(:,2) = [-30;30]; %y coordinates
                               
            %cellposnfile = 'cellposn400cells.mat';
            %cellposnfile = 'cellposnborder.mat';
            cellposnfile = [];
            %cellposnfile = 'headangle_data_H3reverse.mat';
            cellposnfile = 'cellposnOpenCV2H_headangle.mat';


            if exist(cellposnfile,'file') == 2
                disp('Loading cell-position file')
                load(cellposnfile)
            else
                disp('Generating cell-positions')
                cellposn = mbr_cell_distribution(MBRcorners,numcell,celllength);
                
            end
        else
            cellposn = [];
            MBRcorners.cells = [-20 -20;20 20];
            MBRcorners.nocells = [];
        end
        
        htraj = figure;
        hold on
        xlabel('x (um)')
        ylabel('y (um)')
        axis ij
        axis equal
        
        t = CTimeleft(repeatSim);
        for rep = 1:repeatSim
            t.timeleft();
            plotnow = false;
            
            % run simulation
            init = [I A AA AAp YY YYp Mot Run];
            %cellposn = [];
            [timeVec,state] =  MBR_gillespie_func(rxnrate,init,simIterations,attractant,MBRcorners,cellposn);
            
            % compute stats
            if rep == 1
                %simTime = floor(timeVec(end)*0.75);
                simTime = timeVec(end)*0.90;
                fprintf('simTime = %4.4f \n',simTime)
                %simTime = 10; % 10 second simulation time.
            end
            
            [rtRatio,dtheta,dx] = eval_MBR_gillespie(timeVec,state,simTime);
            
            stats.dx(rep,:) = dx;
            stats.dtheta(rep) = dtheta;
            stats.rtratio(rep,:) = rtRatio;
            
            %% Figure of trajectories superimposed %%
            timeIdx = find(timeVec>simTime,1);
            if isempty(timeIdx)
                disp('Sim too short')
                rep
                timeIdx = length(timeVec);
            end
            figure(htraj)
            plot(state.posn(1:timeIdx,1),state.posn(1:timeIdx,2),'-g')
            
        end
        figure(htraj)
        EndTimeIdx = find(timeVec(end)>state.detTime,1,'last') + 1;
        timeaxis = state.detTime(1:EndTimeIdx);
        xDet = state.detPosn(1:EndTimeIdx,1)*10^6;
        yDet = state.detPosn(1:EndTimeIdx,2)*10^6;
        thDet = state.detPosn(1:EndTimeIdx,3);

        hold on
        plot(xDet,yDet,'-b')
        axis ij
        
        %plots
        figure
        hist(stats.dtheta(:))
        %axis([0 100 0 150])
        title('Histogram of angular velocity ratio')
        xlabel('angular velocity')
        ylabel('frequency')
        
        CCWcount = sum(stats.dtheta(:)>0);
        CWcount =  sum(stats.dtheta(:)<0);
        fprintf('CW : CCW = %3.0d : %3.0d \n',CWcount,CCWcount)
        
        meandtheta = mean(stats.dtheta(:));
        sddtheta = std(stats.dtheta(:));
        
        fprintf('Mean : %5.2f  deg/s \n  SD : %5.2f deg/s \n',meandtheta,sddtheta)
        
        
        
        
        %         figure
        %         for i = 1:size(stats.dx,1)
        %             hold on
        %             line([0 stats.dx(i,1)],[0 stats.dx(i,2)])
        %         end
        %         title('Random Walk Planar Displacement')
        %         xlabel('x (um)')
        %         ylabel('y (um)')
        %
        %         avgrtratio = mean(stats.rtratio(:));
        %         stdrtratio = std(stats.rtratio(:));
        %         fprintf('Mean run-tumble ratio: %4.4f, stdev: %4.4f \n',avgrtratio,stdrtratio)
        %
        %         avgdx = mean(stats.dx,1);
        %         stddx = std(stats.dx,1);
        %         fprintf('Mean displacement ratio: (%4.0f,%4.0f), stdev: (%4.0f,%4.0f) \n',avgdx(1),avgdx(2),stddx(1),stddx(2))
        %
    otherwise
        disp('Invalid simMode')
end

%%
% Compute amount of time in Run vs. Tumble
% when Run == 1: run

if simMode<5
    disp('....For the last round')
    [runTime,tumbleTime,rtRatio,dx] = runtumble(timeVec,state,simTime);
    fprintf('Run Time: %4.1f s, Tumble Time: %4.1f s. Run/Tumble: %4.4f \n',runTime, tumbleTime,rtRatio)
    fprintf('Displacement: (%4.0f,%4.0f) um \n',dx(1),dx(2))
    
    tumble = state.chem(1:end-1,8) - state.chem(2:end,8);
    tFreq = sum(tumble==1)/timeVec(end);
    %rFreq = sum(tumble ==-1)/timeVec(end);
    fprintf('Tumble Freq: %4.2f /s \n',tFreq)
end