% aerodonetics.m
%Author: Lenny Martinez

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
%% Reference Solution
r0 = [22; 0; 0; 5]; %[v(0); theta(0); x(0); y(0)]
u = r0;
t = linspace(0,Tmax,N);

options = odeset('Events', @ground_intersection);
[t, u] = ode45(@glider, t, r0, options);


%% Forward Euler Approximation 
dt = .01;
N_fe = t(length(t))/dt;
rfe(1,:) = r0;
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

%% 4-stage Runge-Kutta Approximation
dt_2 = .01;
N_4rk = t(length(t))/dt_2;
rrk(1,:) = r0;
% Start iterative loop
for n = 1:N_4rk
    rk = rrk(n,:);
    trk = n*dt_2;
    %calculate 4 stages
    a = dt_2.*glider(trk, rk);
    b = dt_2.*glider((trk + dt_2/2), (rk + a'./2));
    c = dt_2.*glider((trk + dt_2/2), (rk + b'./2));
    d = dt_2.*glider((n+1)*dt_2, (rk + c'));
    
    % Update using Forward Euler
    rrk(n+1,1) = rrk(n,1) + (1/6)*(a(1) + 2*b(1) + 2*c(1) + d(1));
    rrk(n+1,2) = rrk(n,2) + (1/6)*(a(2) + 2*b(2) + 2*c(2) + d(2));
    rrk(n+1,3) = rrk(n,3) + (1/6)*(a(3) + 2*b(3) + 2*c(3) + d(3));
    rrk(n+1,4) = rrk(n,4) + (1/6)*(a(4) + 2*b(4) + 2*c(4) + d(4));
    
    %Check if glider has hit the ground
    if rrk(n,4) < 0,
        break
    end 
end

%% Plotting Things Part 1: Plot Each Solution Independently

%Plot ODE45 function x-y plane
figure(1); 
plot(u(:,3),u(:,4), 'LineWidth', 2);grid on;
xlabel('x');
ylabel('y');
title('X-Y Trajectory of Lanchester''s Aerodone');

%Plot ODE45 function theta-v plane
figure(2);
plot(u(:,2),u(:,1), 'LineWidth', 2);grid on;
xlabel('Angle of Attack, \theta');
ylabel('Forward velocity, v');
title('Angle of Attack vs. Forward Velocity of Lanchester''s Aerodone');

%Plot Forward Euler Approximation x-y plane
figure(3); 
plot(rfe(:,3),rfe(:,4), 'LineWidth', 2);grid on;
xlabel('x');
ylabel('y');
title('X-Y Trajectory of Forward Euler Approximation of Lanchester''s Aerodone');

%Plot Forward Euler Approximation theta-v plane
figure(4);
plot(rfe(:,2),rfe(:,1), 'LineWidth', 2);grid on;
xlabel('Angle of Attack, \theta');
ylabel('Forward velocity, v');
title('Angle of Attack vs. Forward Velocity of Forward Euler Approximation of Lanchester''s Aerodone');
% 
% %Plot 4-Stage Runge-Kutta Method Approximation x-y plane
% figure(5); 
% plot(rrk(:,3), rrk(:,4) , 'LineWidth', 2); grid on;
% xlabel('x');
% ylabel('y');
% title('X-Y Trajectory of 4-Stage Runge-Kutta Method Approximation of Lanchester''s Aerodone');
% 
% %Plot 4-Stage Runge-Kutta Method Approximation theta-v plane
% figure(6);
% plot(rrk(:,2), rrk(:,1), 'LineWidth', 2); grid on;
% xlabel('Angle of Attack, \theta');
% ylabel('Forward velocity, v');
% title('Angle of Attack vs. Forward Velocity of 4-Stage Runge-Kutta Method Approximation of Lanchester''s Aerodone');

%% Plotting Things Part 2: Comapring Reference Solution to both Forward Euler and RK4 
% 
% Plot x-y plane
figure(7); hold on; grid on 
plot(u(:,3),u(:,4), 'LineWidth', 2) %ODE45 plot
plot(rfe(:,3),rfe(:,4), '--') %Forward Euler approximation
plot(rrk(:,3),rrk(:,4), '*') % 4-Stage Runge-Kutta Method
xlabel('x');
ylabel('y');
title('X-Y Trajectory of Lanchester''s Aerodone')
legend('Reference Solution using ode45 function',['Forward Euler Approximation with \Deltat = ' num2str(dt)],['4-Stage Runge-Kutta Method Approximation with \Deltat = ' num2str(dt_2)], 'LOCATION','BEST') 
% 
% 
% % Plot theta-v plane
% figure(8); hold on;grid on;
% plot(u(:,2),u(:,1), 'LineWidth', 2) %ODE45 plot
% plot(rfe(:,2),rfe(:,1), '--') %Forward Euler approximation
% plot(rrk(:,2),rrk(:,1), '*') % 4-Stage Runge-Kutta Method
% xlabel('Angle of Attack');
% ylabel('Forward velocity (m/s^2)');
% title('Angle of Attack vs. Forward Velocity of Lanchester''s Aerodone')
% legend('Reference Solution using ode45 function',['Forward Euler Approximation with \Deltat = ' num2str(dt)],['4-Stage Runge-Kutta Method Approximation with \Deltat = ' num2str(dt_2)], 'LOCATION','BEST') 


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