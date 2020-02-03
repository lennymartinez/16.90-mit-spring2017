% % nominal values of parameters
clear all
close all
% m       = 6;    % mass of payload and parachute assembly kilograms
% r       = 0.1;  % radius of payload in meters
% Cd      = 0.5;  % payload coefficient of drag
% wx      = 0;    % horizontal wind speed in m/s
% tfree   = 9;    % time before parachute opens
% topen   = 5;    % time it takes for parachute to open completely
% 
% % initial conditions x(0), y(0), vx(0), vy(0)
% u0 = [-340, 500, 50, 0];
% 
% % compute a single trajectory
% [t, u] = payload_sim(u0, m, r, Cd, wx, tfree, topen);
% 
% % plot the trajectory
% close all;
% figure(1);
% plot(u(:,1), u(:,2), '-')
% xlabel('horizontal displacement (m)');
% ylabel('altitude (m)');
% 
% % plot the velocity over time
% figure(2);
% plot(t, sqrt(u(:,3).^2 + u(:,4).^2),'-')
% xlabel('time after release (s)');
% ylabel('velocity (m/s)');

%% Monte Carlo Simulation

% Monte Carlo loop
N = 1000;

% initialize input samples
x     = zeros(N,1);
y     = zeros(N,1);
vx    = zeros(N,1);
vy    = zeros(N,1);
u0    = zeros(N,4);
m     = zeros(N,1);
r     = zeros(N,1);
Cd    = zeros(N,1);
tfree = zeros(N,1);
topen = zeros(N,1);
wx    = zeros(N,1);

% initialize output samples
xf     = zeros(N,1); % landing sites (m)
vf     = zeros(N,1); % impact velocities (m/s)
intact = zeros(N,1); % payload intact after impact (0, 1)


    for i=1:N % consider using parfor, especially for large N
        
        
        % initial conditions
        x(i)     = normrnd(-340, 3);
        y(i)     = normrnd(500, 3);
        vx(i)    = normrnd(50, 0.5);
        vy(i)    = normrnd(0, 0.5);
        
        % random inputs
        m(i)     = trirand(5,6,7);
        r(i)     = trirand(0.09,0.1,0.11);
        Cd(i)    = trirand(0.45,0.5,0.6);
        tfree(i) = trirand(8.95,9,9.05);
        topen(i) = logrand(1.75,0.1);
        wx(i)    = normrnd(0,2);
        
        
        u0(i,:)  = [x(i), y(i), vx(i), vy(i)];
        % simulate trajectory
        [t, u] = payload_sim(u0(i,:), m(i), r(i), Cd(i), wx(i), tfree(i), topen(i));

        % get the impact location and impact velocity
        xf(i) = u(end,1);
        vf(i) = sqrt(u(end,3)^2 + u(end,4)^2);

        % check if payload survived impact
        if rand() < survival(vf(i), 2.7, 0.15),
            intact(i) = 1;
        end

        if i<=20,
            hold on
            plot(u(:,1), u(:,2), '-')
            title('20 Random Trajectories')
            xlabel('Horizontal Displacement (m)');
            ylabel('Altitude (m)');
        end
    end

    % hold off
    % % Plot histograms of input parameters
    % figure()
    % subplot(2,2,1); histogram(x)
    % title('Random x input')
    % xlabel('x (m)')
    % 
    % subplot(2,2,2); histogram(y)
    % title('Random y input')
    % xlabel('y (m)')
    % 
    % subplot(2,2,3); histogram(vx)
    % title('Random v_x input')
    % xlabel('v_x (m/s)')
    % 
    % subplot(2,2,4); histogram(vy)
    % title('Random v_y input')
    % xlabel('v_y (m/s)')
    % 
    % figure()
    % subplot(2,3,1);histogram(m)
    % title('Random Mass input')
    % xlabel('Mass (kg)')
    % 
    % subplot(2,3,2); histogram(r)
    % title('Random radius input')
    % xlabel('radius (m)')
    % 
    % subplot(2,3,3);histogram(Cd)
    % title('Random C_d input')
    % xlabel('C_d')
    % 
    % subplot(2,3,4);histogram(tfree)
    % title('Random t free input')
    % xlabel('t free (s)')
    % 
    % subplot(2,3,5);histogram(topen)
    % title('Random t open input')
    % xlabel('t open (s)')
    % 
    % subplot(2,3,6);histogram(wx)
    % title('Random w_x input')
    % xlabel('w_x (m/s)')
    %% Monte Carlo Analysis

    % Below is starting point for conducting your Monte Carlo analysis

    inside = find(abs(xf) < 1);
    outside = find(abs(xf) > 1);

    survived = find(intact == 1);
    destroyed = find(intact == 0);

    survived_inside = intersect(survived, inside);
    survived_outside = intersect(survived, outside);
    destroyed_inside = intersect(destroyed, inside);
    destroyed_outside = intersect(destroyed, outside);

    % probability of success (outcome A)
    p_A = length(survived_inside) / N;
    p_B = length(survived_outside) / N;
    p_C = length(destroyed_inside) / N;
    p_D = length(destroyed_outside) / N;

    p = [ p_A, p_B; p_C, p_D];
    
% % plot histogram of landing sites of intact payloads
% figure()
% histogram(xf(survived),'BinWidth',5);
% hold on
% histogram(xf(destroyed),'BinWidth',5);
% hold off
% xlabel('x (m)');
% title('Payload Landing Sites')
% legend('Survived', 'Destroyed')
% 
% figure()
% histogram(vf(survived),'BinWidth',1);
% hold on
% histogram(vf(destroyed),'BinWidth',1);
% hold off
% xlabel('v_f (m/s)');
% title('Payload Impact Velocities')
% legend('Survived', 'Destroyed')

%% Importance Sampling Analysis
success = survived_inside;

% f is the indicator function of a successful payload
f = zeros(N,1);
f(success) = 1;

% compute the importance sampling estimate as sum(f.*p./q) / N



