function [dist] = mbr_cell_distribution(corners,numcell,varargin)

%MRB_CELL_DISTRIBUTION randomly determines the position of NUMCELL on and
%mbr of shape defined by corners.
%
% INPUTS: 
% CORNERS: [X0 Y0;X1 Y1; ... ;Xn Yn] (size: nx2)
% Shape of MBR is defined by lines connecting successive corners of the
% MBR, if the first and last point of the input corners are not the same,
% it is assumed that the shape closes by joining the final and first
% coordinate.  (in um)
%
% NUMCELLS: scalar describing how many cells are on the MBR.
%
% VARARGINS: 
% CELLLENGTH: length of cell in um.  Default cell length 5um.
%
% OUTPUT:
% dist: 

% divide MBR into rectangles based on corners

% number of variable inputs
nargin;
if nargin > 0
    celllength = varargin{1};
else
    celllength = 5;
end

%corners = zeros(5,2);
%corners(:,1) = [-20 -20 20 20 -20];
%corners(:,2) = [-20 20 20 -20 -20];

%numcell = 400;

%celllength = 3;

% initialize matrix to store posn of cells
dist = zeros(numcell,3);

% check for corners and make sure shape closes
if corners(end,:) ~= corners(1,:)
    corners = [corners;corners(end,:)];
end

i = 0;
while i < numcell
    intersect = 1;
        while intersect  
            % randomly place a cell
            posnTemp = rand(1,3).*[40,40,360] - [20,20,0];
    
            % check for intersections
            intersect = check_line_intersection(posnTemp,dist,celllength);
            if intersect == 1
                %disp('Regenerate')
                %keyboard
            end
        end
        
    
    % store cell data
    dist(i+1,:) = posnTemp;
    
    % iterate counter
    i = i+1;
    %fprintf('%d cells placed\n',i) 
end

% plot distribution
% determine edge bacteria for a 40um-square
edgecell = zeros(1,numcell); % store 1 if edge bacterium 

% compute location of flagellum
cellangle = dist(:,3);
dbac = celllength*[cosd(cellangle) sind(cellangle)];

bacHead = dist(:,1:2); 
bacTail = bacHead + dbac;

edgecell = max(abs(bacTail),[],2)>20; % setup for 40x40 sq mbr

figure
title('MBR cell distribution')
plot(corners(:,1),corners(:,2),'-k')
for cell = 1:numcell
    hold on
     % determine if it is in the MBR
    hold on
    bacX = [bacHead(cell,1);bacTail(cell,1)];
    bacY = [bacHead(cell,2);bacTail(cell,2)];
    hold on
    if edgecell(cell)
       % OVER EDGE IN RED
       cellplot(cell) = plot(bacX,bacY,'-r','LineWidth',5);
    else
       cellplot(cell) = plot(bacX,bacY,'-b','LineWidth',5);
    end
    
end
axis equal
