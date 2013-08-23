function [edgecell,bacHead,bacTail] = find_edge_bacteria(cellregion,nocellregion,dist,celllength)

% Find Which bacteria are are at the edge.  Returns 1 for edgecells and 0
% for non edge cells.

numcell = size(dist,1);

% compute location of flagellum
cellangle = dist(:,3);
dbac = celllength*[cosd(cellangle) sind(cellangle)];

bacHead = dist(:,1:2); 
bacTail = bacHead + dbac;

% check if bacTail is in a 
% for outside of square
%edgecell = max(abs(bacTail),[],2)>; % setup for 40x40 sq mbr
edgecellregion = is_point_in_box(cellregion,bacTail);

% check if tail is in nocell regions
edgenocellregion = ~is_point_in_box(nocellregion,bacTail);
edgecell = max([edgecellregion,edgenocellregion],[],2);

