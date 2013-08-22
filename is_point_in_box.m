function [notinbox] = is_point_in_box(boxcorners,point)

% Given a box and a set of xy coordinates for point, return a vector inbox
% storing 1 if the point is in the box and 0 if the point is not in the box
%
% INPUT
% define box by bottom left and top right corner SIZE = [2,2] = [x1,y1;x2,y2]
% point SIZE = [n,2] = [x y]
%
% OUTPUT
% inbox SIZE = [n,2]

%% Allow for multiple pairs of corners as well
numpt = length(point);

numbox = size(boxcorners,1)/2;
for i = 1:numbox
    
    inbox = zeros(1,numpt);
    
    %define bounds
    lo = 2*i-1;
    hi= 2*i;
    
    xmin = min(boxcorners(lo:hi,1));
    ymin = min(boxcorners(lo:hi,2));
    xmax = max(boxcorners(lo:hi,1));
    ymax = max(boxcorners(lo:hi,2));
    
    % return 1 if beyond limits
    xyCheck = zeros(numpt,4);
    xyCheck(:,1) = bsxfun(@lt,point(:,1),xmin);
    xyCheck(:,2) = bsxfun(@gt,point(:,1),xmax);
    xyCheck(:,3) = bsxfun(@lt,point(:,2),ymin);
    xyCheck(:,4) = bsxfun(@gt,point(:,2),ymax);
    
    notinbox(:,i) = max(xyCheck,[],2); %search each row and save max value
end

notinbox = max(notinbox,[],2);