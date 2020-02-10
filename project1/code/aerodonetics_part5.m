% aerodonetics_part5.m
% Author: Lenny Martinez
% Using this to answer part 5 of project

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

N = 1001; %number of time steps taken.
Tmax = 100; %max time
t = linspace(0,Tmax,N);


%% Run Forward Euler first to get V, Theta values
dt = .01;
N_fe = t(length(t))/dt;
rfe(1,:) = [29,0,0,10];
% Start iterative loop
for n = 1:N_fe
    %calculate rhs for each variable
    fv = -g*sin(rfe(n,2)) - R_D*(rfe(n,1)^2);
    ftheta = R_L*rfe(n,1) - (g*cos(rfe(n,2)))/rfe(n,1);
    fx = rfe(n,1)*cos(rfe(n,2));
    fy = rfe(n,1)*sin(rfe(n,2));

    % Update using Forward Euler
    rfe(n+1,1) = rfe(n,1) + dt*fv;
    rfe(n+1,2) = rfe(n,2) + dt*ftheta;
    rfe(n+1,3) = rfe(n,3) + dt*fx;
    rfe(n+1,4) = rfe(n,4) + dt*fy;

    if rfe(n,4) < 0,
        break
    end       
end



%% Quiver and Meshgrid
v = linspace(1,30,20);
theta = linspace(-pi,3*pi,20);
[X,Y] = meshgrid(v,theta);
u0 = [0,0];

dv = -g*sin(Y) - R_D*(X.^2);
dtheta = R_L.*X - (g*cos(Y))./X;

% quiver(Y,X,dtheta,dv)
% xlabel('\theta')
% ylabel('v')
% title('Quiver plot for v-\theta')

figure(2)
quiver(Y,X,dtheta,dv); hold on;
xlabel('\theta')
ylabel('v')
title('Quiver plot for v-\theta with solution path for \theta(0)=0, v(0)=0')
plot(0,29, 'r+', 'LineWidth',2)
plot(rfe(:,2),rfe(:,1))
