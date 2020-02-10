% aerodonetics_errors.m
% Author: Lenny Martinez
% Using this to find errors

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

N = 10001; %number of time steps taken.
Tmax = 10; %max time
%% Reference Solution
r0 = [22; 0; 0; 5]; %[v(0); theta(0); x(0); y(0)]
t = linspace(0,Tmax,N);

options = odeset('Events', @ground_intersection,'RelTol', 1e-13, 'AbsTol', 1e-12);
[tref, uref] = ode45(@glider, t, r0, options);


%% Error calculations

dt = [0.05, 0.01, 0.005,0.001];
rfe(1,:) = r0;
rrk(1,:) = r0;
efe = [];
erk = [];

%iterate for different dt values
for j = 1:length(dt),
    N = t(length(t))/dt(j);
    
    % Start iterative loop for Forward Euler
    for n = 1:N
        %calculate rhs for each variable
        fv = -g*sin(rfe(n,2)) - R_D*(rfe(n,1)^2);
        ftheta = R_L*rfe(n,1) - (g*cos(rfe(n,2)))/rfe(n,1);
        fx = rfe(n,1)*cos(rfe(n,2));
        fy = rfe(n,1)*sin(rfe(n,2));

        % Update using Forward Euler
        rfe(n+1,1) = rfe(n,1) + dt(j)*fv;
        rfe(n+1,2) = rfe(n,2) + dt(j)*ftheta;
        rfe(n+1,3) = rfe(n,3) + dt(j)*fx;
        rfe(n+1,4) = rfe(n,4) + dt(j)*fy;

        if rfe(n,4) < 0,
            break
        end       
    end
    
    % Start iterative loop for 4-Stage Runge-Kutta
    for n = 1:N
        rk = rrk(n,:);
        dt_2 = dt(j);
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
    
    %We want error at t = 10
    n10 = 10./dt(j) +1;
    e = max(abs(rfe(n10,:) - uref(end,:)));
    e2 = max(abs(rrk(n10,:) - uref(end,:)));
    efe = [efe, e];
    erk = [erk, e2];
    
end


%% Plotting Errors

% Plot forward euler errors
figure(1)
loglog(dt,efe,'*-'); grid on;
xlabel('Time Step');
ylabel('Global Error');
title('Global Error vs. Time Step for Forward Euler Method');

%plot RK4 errors
figure(2)
loglog(dt,erk,'*-'); grid on;
xlabel('Time Step');
ylabel('Global Error');
title('Global Error vs. Time Step for 4-Stage Runge-Kutta Method');

%compare errors from both methods
figure(3)
loglog(dt,efe,'*-'); grid on; hold on;
loglog(dt,erk,'*-'); grid on;
xlabel('Time Step');
ylabel('Global Error');
legend('Forward Euler Method', '4-Stage Runge-Kutta','LOCATION', 'best');
title('Global Error vs. Time Step ');

%Find Slopes
mdt = log(dt(end)/dt(1));

mfe = log(efe(end)/efe(1))/mdt;
mrk = log(erk(end)/erk(1))/mdt;

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