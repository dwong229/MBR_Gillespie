function PlotSlider_Callback(hObject,event,dxcircle,dycircle,dthcircle,q)

% take slider value for p and update theta plot

slider_value = round(get(hObject,'Value'));

% given that p and q value

% compute dxdt dydt dthdt and update values 
% run deterministic sim
[~,~,~,~,dxdt_body] = runDeterministicModel(1/invkt,1/invkr,p,q,cellposn,edgecell,[0;0;0]');
% returns in m

% unpack velocities
dxdt = dxdt_body(1) * 1e6;
dydt = dxdt_body(2) * 1e6;
dthdt = dxdt_body(3);

% circle [p,q,dxdt]
set(dxcircle,'XData',slider_value,'YData',q,'ZData',dthdt)

% circle [p,q,dydt]


% circle [p,q,dthdt]

