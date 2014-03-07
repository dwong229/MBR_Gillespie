function [dxdt,dydt,dthdt,pplot,qplot] = UpdateVelcityDetModel(pVec,qVec,varargin)

% for changes in kt and kr, update all velocity matrices plotted

% If matrix and plots need to be updated, add handles in varargin:
% order varargin: x, y, dt.

global invkt invkr cellposn edgecell

plotIdx = 0;
for pIdx = 1:length(pVec)
    p = pVec(pIdx);
    for qIdx = 1:length(qVec)
        plotIdx = plotIdx + 1;
        q = qVec(qIdx);
        
        % run deterministic sim
        [~,~,~,~,dxdt_body] = runDeterministicModel(1/invkt,1/invkr,p,q,cellposn,edgecell,[0;0;0]');
        % returns in m
        
        % unpack velocities
        dxdt(plotIdx) = dxdt_body(1) * 1e6;
        dydt(plotIdx) = dxdt_body(2) * 1e6;
        dthdt(plotIdx) = dxdt_body(3);
        
        % assign plot axes
        pplot(plotIdx) = p;
        qplot(plotIdx) = q;
    end
end



if ~isempty(varargin)
    xvelplot = varargin{1};
    yvelplot = varargin{2};
    zvelplot = varargin{3};
    
    set(xvelplot,'ZData',dxdt)
    set(yvelplot,'ZData',dydt)
    set(zvelplot,'ZData',dthdt)
end
    