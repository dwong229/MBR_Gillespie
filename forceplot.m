function [] = forceplot(fighandle,arrowtail,arrowtip)

% input figure handle and plot arrowhead and arrow end to show forces


% make figure handle the current figure
figure(fighandle);
hold on

u = arrowtip(:,1) - arrowtail(:,1);
v = arrowtip(:,2) - arrowtail(:,2);

quiver(arrowtail(:,1),arrowtail(:,2),u,v);
%quiver(X,Y,U,V)