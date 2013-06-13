function [tau,newcellchem,dth] = cell_gillespie(timeVec,rxn_rates,state,varargin)


dth = 0;
if ~isempty(varargin)
    varInput = char(varargin);
end

% unpack rxn
phi = rxn_rates(1);  % I -> A
delta = rxn_rates(2);% A -> I
alpha = rxn_rates(3);% AA -> AAp
beta = rxn_rates(4);% A + AAp + YY -> A + AA + YYp
gamma = rxn_rates(5); % YYp -> YY
mu = rxn_rates(6);   % YYp ->Mot
rho = rxn_rates(7);  % Run

alpha_a = rxn_rates(8);
gamma_a = rxn_rates(9);% A + YYp -> A + YY

% unpack state
I = state(1);
A = state(2);
AA = state(3);
AAp = state(4);
YY = state(5);
YYp = state(6);
Mot = state(7);
Run = state(8);

%% state changes
stateChange = zeros(length(rxn_rates),length(state)); % [rxn rate x state molecules]
stateChange(1,:) = [0 0 1 zeros(1,5)]; % I
stateChange(2,:) = [1 -1 zeros(1,6)]; % A
stateChange(3,:) = [0 0 -1 1 zeros(1,4)]; % AA
stateChange(4,:) = [0 0 1 -1 -1 1 0 0]; % AAp + YY
stateChange(5,:) = [0 0 0 0 1 -1 0 0]; % YYp ->YY
stateChange(6,:) = [zeros(1,6) 0 1]; % YYp -> More likely to tumble, but still running
stateChange(7,:) = [zeros(1,5) 0 0 -1]; % run + YYp -> tumble
stateChange(8,:) = [0 0 -1 1 0 0 0 0]; % CheA -> CheAp
stateChange(9,:) = [0 0 0 0 1 -1 0 0]; % YYp -> YY

rxnRateNow = rxn_rates .* [I A (AA>0) (YY>0 && AAp>0) YYp>0 Run<1 (Run>0 && YYp>0) (A>0 && AAp>0) (A>0 && YYp>0)];

% Compute scaling factor to compute probability of each reaction
rxnScale = 1/sum(rxnRateNow);

rxnProb = rxnScale*rxnRateNow;
rxnBin = cumsum(rxnProb); %generate bound of bins (if  1<r1<2, then rxn = 1)

% Generate random variables: r(1): rxn, r(2): time rxn occurs
r = rand(2,1);

% Compute tau and i

if max(rxnProb)<1
    %disp('More than one rxn possible')
    %       disp('Inequalities')
    probrange(1,:) = [0 rxnBin(1:end-1)];
    probrange(2,:) = rxnBin;
    
    % compute when r(1) is greater then lower bound
    nextRxn = find(r(1)>probrange(1,:),true,'last');
    if r(1)>probrange(2,nextRxn)
        error('probrange error')
    end

else
    %disp('1 rxn possible only')
    nextRxn = find(rxnProb == 1,true,'first');
end
%fprintf('rand(1): %3.3f\n',r(1))
%rxnProb
%rxnBin
%fprintf('Reaction: %3.0f \n',nextRxn)

% Update state
if isequal('MBR',varInput)
    newcellchem = nextRxn;
else
    newcellchem = state + stateChange(nextRxn,:);
end

%%%%%%%%%%%%%%Tumble reorient%%%%%%%%%%%%
if false
% set run state to always be at most 1:
if newcellchem(8) >0
    newcellchem(8) = 1;
elseif newcellchem(8)==0
    if rand(1)>0.5
        dth = normrnd(58,20)*(-1);
    else
        dth = normrnd(58,20)*(1);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute time to event
tau = rxnScale * log(1/r(2));

