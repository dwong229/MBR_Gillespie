function [newstate] = MBRupdatestate(cellstate,nextRxn)

rxn_rate_length = 9;
state_length = 8;

stateChange = zeros(rxn_rate_length,state_length); % [rxn rate x state molecules]
stateChange(1,:) = [0 0 1 zeros(1,5)]; % I
stateChange(2,:) = [1 -1 zeros(1,6)]; % A
stateChange(3,:) = [0 0 -1 1 zeros(1,4)]; % AA
stateChange(4,:) = [0 0 1 -1 -1 1 0 0]; % AAp + YY
stateChange(5,:) = [0 0 0 0 1 -1 0 0]; % YYp ->YY
stateChange(6,:) = [zeros(1,6) 0 1]; % YYp -> More likely to tumble, but still running
stateChange(7,:) = [zeros(1,5) 0 0 -1]; % run + YYp -> tumble
stateChange(8,:) = [0 0 -1 1 0 0 0 0]; % CheA -> CheAp
stateChange(9,:) = [0 0 0 0 1 -1 0 0]; % YYp -> YY

newstate = cellstate + stateChange(nextRxn,:);
