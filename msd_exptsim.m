

transRot = 1; % 1: translation, 2: rotation

if transRot == 1
    disp('Translation')
    load('2Htrackstat.mat')
    xy = H2xy; % nx2 
    th = H2angle; % in deg
else
    disp('Rotation')
    load('H3taackstat.mat')
    xy = H3xy; % nx2
    th = H3angle; % in deg
end

frame2sec = 2; % frames per second

exptTime = frame2sec*[0:length(th)-1];

tau = 2; % amt of time to simulate between updating from expt data

simTime = 0;
% simulate 
while simTime < exptTime(end)  
    
    % update mbr posn and orientation
    
    % simulate for tau time
    mbrState = wrapperMBRgillespiefunc(cellposn,lastState);
     
    % plot update
    
end