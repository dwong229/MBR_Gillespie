function statedot = compute_rxnderivatives(t,state)

% equations setup for incoherent Type 1 AND gate (Sy= 1)
%%%%%%%%%%%%%%%%%%%%
%% parameters of reactions
H = 2;
Betay = 1;
Betaz = 1;
Alphay = 1;
Alphaz = 1;
By = 0.3; %Basal level of Y
Bz = 0;
%Kxz = 1;
%Kxy = 1;
Kyz= 0.5;

Kxy = 1;
Kxz = 1;
%% regulation function
factivator = @(u,K) (u/K)^H/(1+(u/K)^H);
frepressor = @(u,K) 1/(1+(u/K)^H);
% gatefunction transcription factors
% AND-gate
GzAND = @(Xactive,Kxz,Yactive,Kyz) factivator(Xactive,Kxz)*frepressor(Yactive,Kyz);

% OR-gate
% GzOR = 

% protein Z decay rate:
%tstop = 0;
%Zdecay = 

Sy = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xactive = state(1);
Y = state(2);
Z = state(3);
Yactive = Y * Sy;
% time derivatives
Xdot = 0;
Ydot = By + Betay * factivator(Xactive,Kxy) - Alphay*Y;
Zdot = Bz + Betaz * GzAND(Xactive,Kxy,Yactive,Kyz) - Alphaz*Z;

statedot = [Xdot;Ydot;Zdot];