% aerodonetics_trajectories.m
% Author: Lenny Martinez
% Using this to find trajectories for given initial conditions

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

N1 = 401; %number of time steps taken.
Tmax = 100; %max time
t = linspace(0,Tmax,N1);

%All initial condition
u0 = [22,0,0,10];
ua = [29,0,0,10];
ub = [23.1,0,0,10];
uc = [12.0223,-0.0831,0,10];
ud = [6,0,0,10];

dt = 0.01;
N = t(length(t))/dt;

%change this to test the cases!!!!!
rfe(1,:) = ud;

% Start iterative loop for Forward Euler
for n = 1:N
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
    



%% Plotting Trajectories
plot(rfe(:,3),rfe(:,4)); grid on;
xlabel('x');
ylabel('y');
title('X-Y Trajectory for Initial Condition u_d');




%% Functions for Ode45
function [position, isterminal, direction] = ground_intersection(t,u)
    position = u(4); % event triggers when position = 0
    isterminal = 1;  % halt integration
    direction = -1;  % trigger when event function is decreasing
end

function [f] = glider(t, u)
    rho = 1.22; %[kg m^-3]
    g = 9.81; %[m s^-2]
    m = 0.65; %[kg]
    S = 0.06; %[m^2]; Wing Area
    C_D = 0.10; %drag coefficient
    C_L = 1.20; %lift coefficient

    R_D = (rho*C_D*S)/(2*m); %Drag Constant; [m^-1]
    R_L = (rho*C_L*S)/(2*m); %Lift Constant; [m^-1]
    
    
    v = u(1);
    theta = u(2);
    x = u(3);
    y = u(4);
    
    
    dv = -g*sin(theta) - R_D*(v)^2;
    dtheta = R_L*v - ((g*cos(theta))/v);
    dx = v*cos(theta);
    dy = v*sin(theta);
    
    f = [dv; dtheta; dx; dy]; 
end