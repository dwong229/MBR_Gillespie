function QPlotSlider_Callback(hObject,event,dxcircle,dycircle,dthcircle,qtext)

global invkt invkr cellposn edgecell qmin qmax p q

% take slider value for p and update theta plot

slider_value = round(get(hObject,'Value'));
q = slider_value*(qmax-qmin)/100 + qmin;
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
set(dxcircle,'XData',p,'YData',q,'ZData',dxdt)

% circle [p,q,dydt]
set(dycircle,'XData',p,'YData',q,'ZData',dydt)

% circle [p,q,dthdt]
set(dthcircle,'XData',p,'YData',q,'ZData',dthdt)

ppN = p*10^12;
qpN = q*10^12;

fprintf('Forces p:%2.5f pN q:%2.5f pN \n',ppN,qpN)
fprintf('Velocity\nx:%2.3f um/s \ny:%2.3f um/s \nth:%2.3f deg/s \n',dxdt,dydt,dthdt)
disp('--')

% update ptext:
set(qtext,'String',strcat('Select force p: ',num2str(qpN),'pN'))