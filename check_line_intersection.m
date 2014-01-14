function intersect = check_line_intersection(varargin)

%FOR NARGIN = 3
%CHECK_LINE_INTERSECTION(posn,celldist,l) given posn and cell distribution
%and cell length, determines if lines intersect
%
% INPUTS:
% POSN: [X1 Y1 angle(360)]
%
% CELLDIST : (SIZE: n x 3)

%FOR NARGIN = 4
%CHECK_LINE_INTERSECTION([x1a y1a],[x1b y1b],[x2a y2a],[x2b y2b]) given cell posn x1 and y1, check if
%there is an intersection with x2 y2
% 
% INPUTS:
%
% [x1a y1a],[x1b y1b] coordinates of the ends of the cell to check
% intersections, cell_1
%
% [x2a y2a],[x2b y2b] n-rows of cells to check intersections with cell_1.

% OUTPUT:
% Return an array of the index of which cells intersect with cell_1.  

intersect = 0;
troubleshoot = false;

if nargin == 3
    % unpack variables
    posn = varargin{1};
    celldist = varargin{2};
    l = varargin{3};
        
    if sum(celldist(:)) == 0
        % this is the first line placed
        % no intersections
        intersect = 0;
    else
        % compute endpoint of posn
        posnend = posn(1:2) + l*[cosd(posn(3)),sind(posn(3))];
        
        % compute determinants
        detxy = det([posn(1:2);posnend(1:2)]);
        detx = det([posn(1) 1;posnend(1) 1]);
        dety = det([posn(2) 1;posnend(2) 1]);
        
        % number of cells placed
        numcell = find(sum(celldist,2)~=0,1,'last');
        
        % compute end point of cell dist
        cellend = celldist(1:numcell,1:2) + l*[cosd(celldist(1:numcell,3)),sind(celldist(1:numcell,3))];
        
        % find intersections
        i = 1;
        
        if troubleshoot
            plot([posn(1), posnend(1)],[posn(2),posnend(2)],'-r');
            hold on
            h1 = plot([celldist(1,1), cellend(1,1)],[celldist(1,2),cellend(1,2)],'-k');
            axis equal
        end
        
        while intersect == 0 && i <= numcell
            
            if troubleshoot
                set(h1,'XData',[celldist(i,1), cellend(i,1)],'YData',[celldist(i,2),cellend(i,2)]);
            end
            
            detxytemp = det([celldist(i,1:2);cellend(i,1:2)]);
            detxtemp = det([celldist(i,1) 1;cellend(i,1) 1]);
            detytemp = det([celldist(i,2) 1;cellend(i,2) 1]);
            
            % check if cell i intersects
            Px = det([detxy detx;detxytemp detxtemp])/det([detx dety;detxtemp detytemp]);
            Py = det([detxy dety;detxytemp detytemp])/det([detx dety;detxtemp detytemp]);
            
            lim = [min(cellend(i,1),celldist(i,1)),max(cellend(i,1),celldist(i,1))];
            % check if point is on cell segment
            if Px <= lim(2) && Px >= lim(1) && Px <= max(posn(1),posnend(1)) && Px>=min(posn(1),posnend(1))
                % intersection found
                intersect = 1;
                %fprintf('Intersects with cell %d\n',i)
                %keyboard
            end
            i = i+1;
            
        end
        %close all
    end
    
elseif nargin == 4
    % Inputs of endpoint coordinates of cells
    xy1a = varargin{1};
    xy1b = varargin{2};
    xy2a = varargin{3};
    xy2b = varargin{4};
    intersect = [];
%keyboard        
i = 1;
length(xy2b);
    while i <= size(xy2a,1) && length(intersect)<3
        % check for intersection
        cellx = [xy1a(1) xy1b(1)];
        celly = [xy1a(2) xy1b(2)];
        checkx = [xy2a(i,1) xy2b(i,1)];
        checky = [xy2a(i,2) xy2b(i,2)];
        
        [Px,Py] = polyxpoly(cellx,celly,checkx,checky);
        
        % limits defined at [minx miny maxx maxy] for each cell
        celllim = [min(cellx),min(celly),max(cellx),max(celly)];
        checklim = [min(checkx),max(checky),max(checkx),max(checky)];
        
        % check if point is on cell segment
        if ~isempty(Px)
            disp('Intersection detected')
            
            if Px <= celllim(3) && Px >= celllim(1) && Px <= checklim(3) && Px >= checklim(1)
                % intersection found
                %keyboard
                intersect = [intersect i];
                %fprintf('Intersects with cell %d\n',i)
                
                %keyboard
            end
        end
        i = i+1;
    end
    
else
    disp('function: CHECK_LINE_INTERSECTION - Invalid number of inputs')
end