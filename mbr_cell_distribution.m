function [dist,edgecell,bacHead,bacTail] = mbr_cell_distribution(corners,numcell,varargin)

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

% VARARGIN
% {1} : celllength
% {2} : plottroubleshoot

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
%if corners.cells(end,:) ~= corners.cells(1,:)
%    corners.cells = [corners.cells;corners.cells(end,:)];
%end

xrange = range(corners.cells(:,1));
yrange = range(corners.cells(:,2));

i = 0;
while i < numcell
    intersect = 1;
    while intersect
        % randomly place a cell
        posnTemp = rand(1,3).*[xrange,yrange,360] + [corners.cells(1,1),corners.cells(1,2),0];
        
        % check if cell positions is valid (e.g. not in no cell region square)
        if ~isempty(corners.nocells)
            if is_point_in_box(corners.nocells,posnTemp(1,1:2))==0
                % point in no cell region
                intersect = 1;
                %                disp('point in nocell region')
                %                keyboard
            else
                % point not in no cell region, check for intersections
                intersect = check_line_intersection(posnTemp,dist,celllength);
            end
        else
            % no need "nocell" region to check, just
            % check for intersections
            intersect = check_line_intersection(posnTemp,dist,celllength);
            %if intersect == 1
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
%% EDGE BACTERIA
% determine edge bacteria
edgecell = zeros(1,numcell); % store 1 if edge bacterium

[edgecell,bacHead,bacTail] = find_edge_bacteria(corners.cells,corners.nocells,dist,celllength);


%% PLOT
%plottroubleshoot = false;
if nargin == 4 
   plottroubleshoot = varargin{2};
else
    plottroubleshoot = false;
end

if plottroubleshoot
figure
title('MBR cell distribution')
x1 = corners.cells(1,1);
y1 = corners.cells(1,2);
x2 = corners.cells(2,1);
y2 = corners.cells(2,2);

cornerallpts(:,1) = [x1,x2,x2,x1,x1];
cornerallpts(:,2) = [y1,y1,y2,y2,y1];

plot(cornerallpts(:,1),cornerallpts(:,2),'-k')

if ~isempty(corners.nocells)
    % plot no cell region
    x1 = corners.nocells(1,1);
    y1 = corners.nocells(1,2);
    x2 = corners.nocells(2,1);
    y2 = corners.nocells(2,2);
    
    hold on
    cornerallpts(:,1) = [x1,x2,x2,x1,x1];
    cornerallpts(:,2) = [y1,y1,y2,y2,y1];
    plot(cornerallpts(:,1),cornerallpts(:,2),'-k')
end

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
end
axis equal
