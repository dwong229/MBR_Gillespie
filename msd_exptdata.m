%% Compute MSD for tracked MBRs

load('2Htrackstat.mat')
% H2angle (in deg)
% H2xy (in microns)

x = H2xy(:,1);
y = H2xy(:,2);


d = 1;

i = 1;
notatend = true;

while notatend
    msd{i} = sum((H2xy(end-i+1:-1:1+i,:) - H2xy(end-1:-1:1,:)).^2,2);
    
    
    notatend = false;
end