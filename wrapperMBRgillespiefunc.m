function [timeVec,state] = wrapperMBRgillespiefunc(cellposn,MBRcorners,lastState,tau)

% INPUT: 
% 
% OUTPUT: 
% 

% Simulation parameters
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

I = 1; % no chemoattractant
A = 0;
AA = 0;  %CheA
AAp = 0; %CheA-p
YY = 1;  %CheY
YYp = 1; %CheY-p
Mot = 1; %CheY-p attached to motor
Run = 1;

init.chem = [I A AA AAp YY YYp Mot Run];

% initialize posn
init.dyn = lastState; %start at [x,y,theta]

attractant = @(x) 0;

simIterations = 1000 * (tau + 1);

simTimeLongEnough = true;

while simTimeLongEnough
    
    %% Run MBR gillespie func
    [timeVec,state] =  MBR_gillespie_func(rxnrate,init,simIterations,attractant,MBRcorners,cellposn);
    disp('End Time')
    timeVec(end)
    
    if timeVec(end) < tau
        simIterations = ceil(simIterations * tau/timeVec(end));
    else
        simTimeLongEnough = false;
    end
    
end

        
