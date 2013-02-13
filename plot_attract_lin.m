function [] = plot_attract_lin(fighandle,path)

figure = fighandle;
% figure out x/y axis bounds:
hold on

[~,xs,ys] = attract_lin([0,0]);
xmin = path(1)-10;
xmax = path(2)+10;
ymin = path(3)-10;
ymax = path(4)+10;
%xcoord = linspace(min(xs,path(1))-10,max(xs,path(2))+10,50);
%ycoord = linspace(min(ys,path(3))-10,max(ys,path(4))+10,50);
xcoord = linspace(path(1)-10,path(2)+10,50);
ycoord = linspace(path(3)-10,path(4)+10,50);
% copied from attract_exp.m
yconc = 1 - abs(ys - ycoord)./(ys*2);

concgrad = repmat(yconc',[1,length(xcoord)])

%contourf(xcoord,ycoord,concgrad)
pcolor(xcoord,ycoord,concgrad)
axis([xmin, xmax,ymin,ymax])
%mesh(xcoord,ycoord,concgrad)        

%concgrad = arrayfun(@(xcoord,ycoord)(exp(-sqrt((xs-xcoord)^2+(ys-ycoord)^2))),xcoord,ycoord);
shading flat
shading interp
