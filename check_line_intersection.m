function intersect = check_line_intersection(posn,celldist,l)

%CHECK_LINE_INTERSECTION(plot,celldist,l) given posn and cell distribution
%and cell length, determines if lines intersect
%
% INPUTS: 
% POSN: [X1 Y1 angle(360)]
% 
% CELLDIST : (SIZE: n x 3)
intersect = 0;

% 
troubleshoot = false;

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