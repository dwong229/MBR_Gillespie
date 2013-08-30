idx = -19:2:19;
numcell = length(idx);
cellposn = zeros(numcell*4,3);

% y = 20, deg = 90
cellposn(1:1*numcell,:) = [idx',ones(numcell,1)*20,ones(numcell,1)*90];

% y = -20, deg = 270
cellposn(1*numcell + 1:2*numcell,:) = [idx',ones(numcell,1)*(-20),ones(numcell,1)*270];

% x = -20, deg = 180;
cellposn(2*numcell + 1:3*numcell,:) = [ones(numcell,1)*(-20),idx',ones(numcell,1)*180];


% x = 20, deg = 0
cellposn(3*numcell + 1:4*numcell,:) = [ones(numcell,1)*(20),idx',ones(numcell,1)*0];

save('cellposnborder.mat','cellposn')