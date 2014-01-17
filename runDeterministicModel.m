function [MBRx,MBRy,MBRth] = runDeterministicModel(kt,kr,runforce,sideforce,cellposn)

% inputs

% outputs


%Simulation of Serratia Marcescens blotted on SU8 rectangular substrate
%Based on Modeling, control and experimental characterization of
%microbiorobots - Selman Sakar et al. (2011)

%System : 50um x 100um, B1 = 13.03um/(s pN), B2 = -43.64um/(s pN), 
%B3 = 1.24 (rad/s pN)

close all
clear
clc

time = 65000;
timestep = 1/1000;

%distances in micro meters (um), forces in pic Newtons (pN)
%Define Constants

width = 50; %um
length = 100; %um

B1 = 1/kt ;%um/(s pN)

B2 = -15.64; %um/(s pN)
B3 = 0.24; % rad/(s pN) 

rho = 0:.1:1.2; %cont * Nb,e/Nb
%rho = 0.1:.1:1.1;
%rho = [0:0.1:1 1.05 1.1 1.15];

% constant from tangential force
%B4edge = -5; % rad/(s pN)

%transition rate
ts(1) = 1; %1/s
ts(2) = 10; %10/s
%force exerted at different states [N]
p(1) = runforce; %e-9 %run  [pN]
p(2) = 0; %tumble
%number of bacterium
n = size(cellposn,1);

%expectation of force exerted
pbar = ts(2)/(ts(1)+ts(2))*p(1)+ts(1)/(ts(1)+ts(2))*p(2)

qbar = sideforce; %(pN)

G1 = rho(k)*8;
G2 = rho(k)*8;
G3 = rho(k)*(-1.7);

%Initialize posn
r1 = [0;0;0]; %r(1) = r_x, r(2)=r_y, r(3) = phi
r2 = [0;0;0];
rhistory = zeros(3,time);

for i = 1:time
xdot = (pbar *B1+qbar*G1)*cos(r1(3)) - (pbar*B2+qbar*G2)*sin(r1(3));
ydot = (pbar *B1+qbar*G1)*sin(r1(3)) + (pbar*B2+qbar*G2)*cos(r1(3));
phidot = pbar*B3 + qbar*G3;

r2(1) = r1(1) + (timestep)*xdot;
r2(2) = r1(2) + (timestep)*ydot;
r2(3) = r1(3) + (timestep)*phidot;

rhistory(:,i) = r2;
r1 = r2;

end

timeaxis = timestep*[1:time];
stat = strcat('(',num2str(rho(k)),',',num2str(phidot,'%3.2f'),')');

hold on
plot(rhistory(1,1:60000),rhistory(2,1:60000),'-b','LineWidth',3)
plot(rhistory(1,1),rhistory(2,1),'go',rhistory(1,60000),rhistory(2,60000),'rx','MarkerSize',10)    
%text(rhistory(1,end),rhistory(2,end),stat,'HorizontalAlignment','center','BackgroundColor',[.7 .9 .7],'FontSize',15)

xlabel('x [\mum]','fontsize',20)
ylabel('y [\mum]','fontsize',20)
axis equal

figure
subplot(3,1,1)
plot(timeaxis,rhistory(1,:))
ylabel('x','fontsize',20)

subplot(3,1,2)
plot(timeaxis,rhistory(2,:))
ylabel('y','fontsize',20)

subplot(3,1,3)
plot(timeaxis,rhistory(3,:))
ylabel('\phi','fontsize',20)
xlabel('Time [s]','fontsize',20)

G3 = -1.7;
figure
phidot = pbar*B3 + qbar*G3*rho;
plot(rho,phidot,'-sb')
xlabel('\rho','fontsize',20)
ylabel('$\dot{\phi}$', 'interpreter', 'latex','fontsize',20)
%----------------------------------------------