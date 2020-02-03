% aerodonetics_part6.m
% Author: Lenny Martinez
% Using this to work on part 6

close all; clear all;

%% Defining Constants and parameters
rho = 1.22; %[kg m^-3]
g = 9.81; %[m s^-2]
m = 0.65; %[kg]
S = 0.06; %[m^2]; Wing Area
C_D = 0.10; %drag coefficient
C_L = 1.20; %lift coefficient

R_D = (rho*C_D*S)/(2*m); %Drag Constant; [m^-1]
R_L = (rho*C_L*S)/(2*m); %Lift Constant; [m^-1]

%% Quiver and Meshgrid
v = linspace(1,30,50);
theta = linspace(-pi,3*pi,50);
[X,Y] = meshgrid(v,theta);
dv = -g*sin(Y) - R_D*(X.^2);
dtheta = R_L.*X - (g*cos(Y))./X;


%% Part 6
th0 = atan(-R_D/R_L);
v0 = sqrt(cos(th0)/R_L);
f0 = [-g*sin(th0) - R_D*(v0^2); R_L*v0 - (g*cos(th0))/v0];
dt = .01;
rfe(1,:) = [v0,th0,0,0.75];
% Start iterative loop
for n = 1:100
    %calculate rhs for each variable
    vlin = f0(1) -g*cos(th0)*rfe(n,2) -2*R_D*v0*rfe(n,1);
    thlin = f0(2) + ((g*sin(th0))/v0)*rfe(n,2) + (R_L +g*(cos(th0)/(v0^2)))*rfe(n,1);
    ftheta = R_L*rfe(n,1) - (g*cos(rfe(n,2)))/rfe(n,1);
    fx = rfe(n,1)*cos(rfe(n,2));
    fy = rfe(n,1)*sin(rfe(n,2));
    
    % Update using Forward Euler
    rfe(n+1,1) = rfe(n,1) + dt*vlin;
    rfe(n+1,2) = rfe(n,2) + dt*thlin;
    rfe(n+1,3) = rfe(n,3) + dt*fx;
    rfe(n+1,4) = rfe(n,4) + dt*fy;

    if rfe(n,4) < 0,
        break
    end       
end

quiver(Y,X,dtheta,dv), hold on; grid on;
plot(th0,v0,'*','LineWidth', 3);
plot(0,0,'+')
xlabel('\theta')
ylabel('v')
title('Quiver plot with \theta^*, v^* plotted')
