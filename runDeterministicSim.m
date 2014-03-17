function [dx dy dth] = runDeterministicSimVelocity(celldist,p,q)

% function that takes cell distribution,p,q, and returns velocity in body frame 
% INPUT
% celldist: contains info on cell position and edge cells
% p: propulsion force
% q: side wall force

% OUTPUT
% dx: velocity in the x direction in um/s
% dy: velocity in the y direction in um/s
% dth: velcity in the theta direction in deg/s

