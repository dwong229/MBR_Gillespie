function dragSlider_Callback(hObject,event,rOrT,texthandle,xvelplot,yvelplot,thvelplot,dxcircle,dycircle,dthcircle)

global invkt invkr cellposn edgecell krmin krmax ktmin ktmax p q pVec qVec

% for drag slider update.  

slider_value = round(get(hObject,'Value'));

disp(slider_value)

if strcmp(rOrT,'r')
    % update rotation drag
    %disp('updating kr')
    kr = slider_value*(krmax-krmin)/100 + krmin;
    invkr = 1/kr;
    set(texthandle,'String',strcat('Select drag inv(kr): ',num2str(invkr,4),'m/Ns'))
else strcmp(rOrT,'t')
    %disp('updating kt')
    % update translation drag
    kt = slider_value*(ktmax-ktmin)/100 + ktmin;
    invkt = 1/kt;
    set(texthandle,'String',strcat('Select drag inv(kt): ',num2str(invkt,4),'m/Ns'))
end



% Update plots
[dxdt,dydt,dthdt,pplot,qplot] = UpdateVelcityDetModel(pVec,qVec);

% update _velplot

set(xvelplot,'XData',pplot,'YData',qplot,'ZData',dxdt)
set(yvelplot,'XData',pplot,'YData',qplot,'ZData',dydt)
set(thvelplot,'XData',pplot,'YData',qplot,'ZData',dthdt)

% update position of x:
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

%% plot display
ppN = p*10^12;
qpN = q*10^12;

fprintf('Forces p:%2.5f pN q:%2.5f pN \n',ppN,qpN)
fprintf('Velocity\nx:%2.3f um/s \ny:%2.3f um/s \nth:%2.3f deg/s \n',dxdt,dydt,dthdt)
disp('--')

%% Update mbr Vector plot:
updateVector(dxdt,dydt,dthdt)
