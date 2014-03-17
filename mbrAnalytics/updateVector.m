function [] = updateVector(dx,dy,dth)

% given update of velocities, update mbr velocity vector

global vvector thvector

% normalize
%normv = [dx dy]/norm([dx dy])*15
normv = [dx dy]/10;

set(vvector,'UData',normv(1),'VData',normv(2))
set(thvector,'VData',dth)