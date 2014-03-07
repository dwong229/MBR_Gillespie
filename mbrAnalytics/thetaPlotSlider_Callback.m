function thetaPlotSlider_Callback(hObject,event,pplot,qplot,varargin)

% take slider value for p and update htheta plot location of circle
% INPUT
% hObject: plot3(x,y,dthdt), circle around point (x,y,dthdt)
% pplot,qplot,dthdt 
% 


slider_value = round(get(hObject,'Value'));

% given that p and q value
 