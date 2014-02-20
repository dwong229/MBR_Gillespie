% function to determine if cell is part of MBR or not



% INPUT
%   - cell distribution [x1 y1 x2 y2]
%   - corners of MBR


% OUTPUT
%   - cellposn : [x y th] posn of head (x,y) 
%   

clear all
close all

%% Importing data
%cellposnfilename = 'cellposnOpenCV2H.txt'; % Translations 
datafile = '2H';
cellposnfilename = 'cellposnOpenCVH3c.txt'; % Rotation 
datafile = 'H3';

fid = fopen(cellposnfilename,'r');

% define corners of MBR:
rawdata = fscanf(fid,'%f %f %f %f');
fclose(fid);

cellposnends = zeros(length(rawdata)/4,4)';
cellposnends(:) = rawdata;
cellposnends = cellposnends';

%% Troubleshoot: plot lines
% figure()
% hold all
% for i = 2:length(cellposnends)
%     line(cellposnends(i,[1,3]),cellposnends(i-1,[2,4]));
% end
% 
% axis ij
% axis image

%% Transform coordinates into 

if datafile == '2H'
    % Define coordinates of corners of H from pixels to um and rotated to
    % aligned (from headstail2headangle.m:
    MBRcenter = [272.7,270];
    
    LLCorner = [40 62];
    LRCorner = [471 33];
    %MBRangle = rad2deg(atan2(33-62,471-40))
    MBRangle = -2;
    pix2mic = 494/60; % pixels per 60um
    
    % just rotating without moving coordinates is not quite correct...
    rotMat = @(x)[cosd(x) -sind(x);sind(x) cosd(x)];
    TR = zeros(3,3);
    TR(1:2,1:2) = rotMat(-MBRangle);
    MBRcenter = rotMat(-MBRangle)*MBRcenter';
    TR(:,3) = [-MBRcenter;1];
    
    nCells = length(cellposnends)*2;
    allcellcoords = [cellposnends(:,1:2);cellposnends(:,3:4)];
    homoCoord = [allcellcoords';ones(1,nCells)];
    homoStraight = [TR*homoCoord]';
    
    % convert pixels to microns
    allcellcoord_data = homoStraight(:,1:2)/pix2mic;
    
    cellposnends_microns = [allcellcoord_data(1:nCells/2,1:2) allcellcoord_data(nCells/2+1:nCells,1:2)]
    
    % Troubleshoot figure
    % figure()
    % hold all
    % for i = 2:length(cellposnends)
    %     line(cellposnends_microns(i,[1,3]),cellposnends_microns(i,[2,4]));
    % end
    %
    % axis ij
    % axis image
    
    %% Determine which cells are within MBR
    % check 3 different rectangles for 40um H
    
    % corner pixels from image
    cornersPixels =    [40.3234  519.0640;
        16.5570   62.0286;
        128.7755   52.8470;
        151.2192  205.8723;
        385.8579  194.6505;
        371.5756   39.5849;
        485.8345   24.2823;
        524.6009  484.3783;
        415.4428  489.4792;
        396.0596  324.2119;
        139 335;
        175  505];
    % 11: 160.4007  342.5749;
    
elseif datafile == 'H3'
    % Define coordinates of corners of H from pixels to um and rotated to
    % aligned (from headstail2headangle.m:
    MBRcenter = [272.7,270];
    
    LLCorner = [40 62];
    LRCorner = [471 33];
    %MBRangle = rad2deg(atan2(33-62,471-40))
    MBRangle = -30;
    pix2mic = 494/60; % pixels per 60um
    
    % just rotating without moving coordinates is not quite correct...
    rotMat = @(x)[cosd(x) -sind(x);sind(x) cosd(x)];
    TR = zeros(3,3);
    TR(1:2,1:2) = rotMat(-MBRangle);
    MBRcenter = rotMat(-MBRangle)*MBRcenter';
    TR(:,3) = [-MBRcenter;1];
    
    nCells = length(cellposnends)*2;
    allcellcoords = [cellposnends(:,1:2);cellposnends(:,3:4)];
    homoCoord = [allcellcoords';ones(1,nCells)];
    homoStraight = [TR*homoCoord]';
    
    % convert pixels to microns
    allcellcoord_data = homoStraight(:,1:2)/pix2mic;
    
    cellposnends_microns = [allcellcoord_data(1:nCells/2,1:2) allcellcoord_data(nCells/2+1:nCells,1:2)]
    
    % Troubleshoot figure
    % figure()
    % hold all
    % for i = 2:length(cellposnends)
    %     line(cellposnends_microns(i,[1,3]),cellposnends_microns(i,[2,4]));
    % end
    %
    % axis ij
    % axis image
    
    %% Determine which cells are within MBR
    % check 3 different rectangles for 40um H
    
    % corner pixels from image
    cornersPixels =    [40.3234  519.0640;
        16.5570   62.0286;
        128.7755   52.8470;
        151.2192  205.8723;
        385.8579  194.6505;
        371.5756   39.5849;
        485.8345   24.2823;
        524.6009  484.3783;
        415.4428  489.4792;
        396.0596  324.2119;
        139 335;
        175  505];
    
    
    
    
    
end

homoCoord = [cornersPixels';ones(1,12)];
homoStraight = [TR*homoCoord]';

% convert pixels to microns
corners = homoStraight(:,1:2)/pix2mic;

rect1 = [corners(1,:);corners(3,:)];
rect2 = [corners(11,:);corners(5,:)];
rect3 = [corners(8,:);corners(6,:)];
figure
hold on
pos_data = [];
j = 1;
for i = 1:length(cellposnends_microns)
    if is_point_in_box(rect1,cellposnends_microns(i,1:2))==0
        %line(cellposnends_microns(i,[1,3]),cellposnends_microns(i,[2,4]));
        pos_data(j,:) = cellposnends_microns(i,:);
        j = j+1;
    elseif is_point_in_box(rect2,cellposnends_microns(i,1:2))==0
        %line(cellposnends_microns(i,[1,3]),cellposnends_microns(i,[2,4]));
        pos_data(j,:) = cellposnends_microns(i,:);
        j = j+1;
    elseif is_point_in_box(rect3,cellposnends_microns(i,1:2))==0
        %line(cellposnends_microns(i,[1,3]),cellposnends_microns(i,[2,4]));
        pos_data(j,:) = cellposnends_microns(i,:);
        j = j+1;
    end   
end
axis equal 
axis ij

figure
xCellAll = pos_data(:,[1,3]);    
yCellAll = pos_data(:,[2,4]);
l2 = line(xCellAll',yCellAll')
set(l2,'color','b');

%% post process selection such that if cell overlap, average them into 2 cell:
% use pos_data to keep track of which cells have not been checked
% use filtered_pos_data to keep track of which cells have not been checked
i = 1;
j = 1; %idx for new averaged coord counter
filtered_pos_data = [];
pos_data(:,5) = [1:length(pos_data)]; % keep track of l2 idx
% pos_data = pos_data(1:5,:)
% 
%     xCellTemp = [pos_data(1,[1,3])];
%     yCellTemp = [pos_data(1,[2,4])];
%     
%     l1 = line(xCellTemp,yCellTemp);
%     set(l1,'Color','r')
%     
%     xCellAll = pos_data(2:end,[1,3]);    
%     yCellAll = pos_data(2:end,[2,4]);
%     
%     l2 = line(xCellAll',yCellAll');
%     set(l2,'color','b');
%     

% set line color to be red if replaced
%                      green if deleted
%                       blue if no change

while size(pos_data,1)>1
%     
     xCellTemp = [pos_data(1,[1,3])];
     yCellTemp = [pos_data(1,[2,4])];
%     
%     %set(l1,'XData',xCellTemp,'YData',yCellTemp,'Color','r');
%     l1 = line(xCellTemp,yCellTemp);
%     set(l1,'Color','r')
%     
     xCellAll = pos_data(2:end,[1,3]);    
     yCellAll = pos_data(2:end,[2,4]);
%     
%     l2 = line(xCellAll',yCellAll');
%     set(l2,'color','b');
%     
%     drawnow
%     
    % check for number of intersections   
    intersection = check_line_intersection(pos_data(1,1:2),pos_data(1,3:4),pos_data(2:end,1:2),pos_data(2:end,3:4));
    
    if isempty(intersection)
        %disp('no intersection')
        filtered_pos_data = [filtered_pos_data;pos_data(1,:)];      
        
    elseif length(intersection) ==1
        % if cell intersects with one other cell, average the values for length and angle
        %disp('2 cells intersect')
        
        coords = pos_data(intersection+1,:);
        distmat = dist([pos_data(1,[1:2])',coords(1:2)' coords(3:4)']);
        set(l2(pos_data(1,5)),'Color','y')
        set(l2(pos_data(intersection+1,5)),'Color','y')
        
        idx = 1;
        if distmat(1,3) < distmat(1,2)
            idx = 2;
            % coordinates for pt1
            %xCellTemp
            %yCellTemp
            %coords
            newcoords = [mean([xCellTemp(1),coords(3)]),mean([yCellTemp(1),coords(4)]), mean([xCellTemp(2),coords(1)]), mean([yCellTemp(2),coords(2)]) 0];
            % coordinates for pt2
        else
            %xCellTemp
            %yCellTemp
            %coords
            newcoords = [mean([xCellTemp(1),coords(1)]), mean([yCellTemp(1),coords(2)]), mean([xCellTemp(2),coords(3)]), mean([yCellTemp(2),coords(4)]) 0];
        end
        l3(j) = line(newcoords([1,3]),newcoords([2,4]));
        set(l3(j),'Color','r')
        j = j+1;
        filtered_pos_data = [filtered_pos_data;newcoords];
                
        % delete row:
        pos_data(intersection + 1,:) = [];
    else
        % if cell interects with more than 1 other cell throw it away
        %disp('More than 2 cells intersect')
        set(l2(pos_data(1,5)),'Color','g')
        
    end
    
    pos_data(1,:) = [];
    
end
if size(pos_data,1) > 0
    filtered_pos_data = [filtered_pos_data;pos_data(1,:)];
end
axis equal 
axis ij
keyboard
pos_data = filtered_pos_data(:,1:4);
% Determine posn of head
%% for translation
% Compute angle of bacteria
% return cellposn in headangle matrix
for i = 1:length(pos_data)
    
    % store temp variables for pairs of coordinates for cells
    x = pos_data(i,[1,3]);
    y = pos_data(i,[2,4]);
    
    cellLength(i) = pdist([x' y']);
    
    % head is coord where x is more negative
    if x(1)<x(2)
        headIdx = 1;
        disp('case 1')
        th = rad2deg(atan2(diff(y),diff(x)));
    else
        headIdx = 2;
        disp('case 2')
        % add - to diff's to flip head and tail
        th = rad2deg(atan2(-diff(y),-diff(x))); 
    end
    
    % save coordinate of head and th in degrees
    headangle_data(i,:) = [x(headIdx),y(headIdx),th];
    line(x,y);
    plot(x(headIdx),y(headIdx),'.r')
end

cellposn = headangle_data;
save( 'cellposnOpenCVH3_headangle.mat','cellposn')
