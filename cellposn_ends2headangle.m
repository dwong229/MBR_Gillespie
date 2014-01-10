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
cellposnfilename = 'cellposnOpenCV2H.txt';
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
% check 3 different rectangles
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
homoCoord = [cornersPixels';ones(1,12)];
homoStraight = [TR*homoCoord]';

% convert pixels to microns
corners = homoStraight(:,1:2)/pix2mic;

rect1 = [corners(1,:);corners(3,:)];
rect2 = [corners(11,:);corners(5,:)];
rect3 = [corners(8,:);corners(6,:)];
figure
hold on
for i = 1:length(cellposnends_microns)
    if is_point_in_box(rect1,cellposnends_microns(i,1:2))==0
        line(cellposnends_microns(i,[1,3]),cellposnends_microns(i,[2,4]));
    elseif is_point_in_box(rect2,cellposnends_microns(i,1:2))==0
        line(cellposnends_microns(i,[1,3]),cellposnends_microns(i,[2,4]));
    elseif is_point_in_box(rect3,cellposnends_microns(i,1:2))==0
        line(cellposnends_microns(i,[1,3]),cellposnends_microns(i,[2,4]));
    end   
end
axis equal 
axis ij

% Determine posn of head


% Compute angle of bacteria


% return cellposn in headangle matrix
