function thetaPlotSlider_Callback(hObject,event,pplot,qplot,varargin)

% take slider value for p and update theta plot

slider_value = round(get(hObject,'Value'));

% given that p and q value
 