%% Monte Carlo Simulation

close all; clear all;

% Monte Carlo loop
N = 1000;
load 'biasing_dist'
% draw N samples from the biasing distribution as a Nx10 matrix
samples = mvnrnd(biasing_dist.mean, biasing_dist.cov, N);

% compute pdf of the biasing distribution q(...) so that q is a Nx1 vector 
% where q(i) is the density evaluated at samples(i,:)
q = mvnpdf(samples, biasing_dist.mean, biasing_dist.cov);

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
p     = zeros(N,10);
u0    = zeros(N,4);
% initialize output samples
xf     = zeros(N,1); % landing sites (m)
vf     = zeros(N,1); % impact velocities (m/s)
intact = zeros(N,1); % payload intact after impact (0, 1)


    for i=1:N % consider using parfor, especially for large N

        % Draw samples from the input distributions
        % unpack row i of the samples

        m(i)       = samples(i,1);
        r(i)       = samples(i,2);
        Cd(i)      = samples(i,3);
        wx(i)      = samples(i,4);
        tfree(i)   = samples(i,5);
        topen(i)   = samples(i,6); 
        x(i)       = samples(i,7);
        y(i)       = samples(i,8);
        vx(i)      = samples(i,9);
        vy(i)      = samples(i,10);
        
        
        % TODO: evaluate pdf of "nominal" distribution p(...) at samples(i,:)
        % Evaluate all log-normal random variables: topen
        p(i,5) = (1/(topen(i)*0.15*sqrt(2*pi())))*exp((log(topen(i))-1.65)^2/(-2*(0.15^2)));
        
        % Evaluate all the normal random variables: x,y,vx,vy,wx
        p(i,7) = (1/sqrt(2*pi()*(3.0^2)))*exp((x(i)+340)^2/(-2*(3.0^2)));
        p(i,8) = (1/sqrt(2*pi()*(3.0^2)))*exp((y(i)-500)^2/(-2*(3.0^2)));
        p(i,9) = (1/sqrt(2*pi()*(0.5^2)))*exp((vx(i)-50)^2/(-2*(0.5^2)));
        p(i,10) = (1/sqrt(2*pi()*(0.5^2)))*exp((vy(i)-0)^2/(-2*(0.5^2)));
        p(i,4) = (1/sqrt(2*pi()*(2.0^2)))*exp((wx(i)-0)^2/(-2*(2.0^2)));
        
        %Evaluate all the triangular random variables: m,r,Cd,tfree
        j = m(i); a= 4.0; b = 6.0; c= 9.0;
        if x < a
            p(i,1) = 0;
        elseif x > c
            p(i,1) = 0;
        elseif x >= a && x < b
            p(i,1) = (2*(x-a))/((c-a)*(b-a));
        elseif x == b
            p(i,1) = 2/(c-a);
        else
            p(i,1) = 2*(c-x)/((c-a)*(c-b));
        end
        
        j = r(i); a= 0.09; b = 0.10; c= 0.11;
        if x < a
            p(i,2) = 0;
        elseif x > c
            p(i,2) = 0;
        elseif x >= a && x < b
            p(i,2) = (2*(x-a))/((c-a)*(b-a));
        elseif x == b
            p(i,2) = 2/(c-a);
        else
            p(i,2) = 2*(c-x)/((c-a)*(c-b));
        end
                
        j = Cd(i); a= 0.45; b = 0.5; c= 0.60;
        if x < a
            p(i,3) = 0;
        elseif x > c
            p(i,3) = 0;
        elseif x >= a && x < b
            p(i,3) = (2*(x-a))/((c-a)*(b-a));
        elseif x == b
            p(i,3) = 2/(c-a);
        else
            p(i,3) = 2*(c-x)/((c-a)*(c-b));
        end
            
        j = tfree(i); a= 8.95; b = 9; c= 9.05;
        if x < a
            p(i,5) = 0;
        elseif x > c
            p(i,5) = 0;
        elseif x >= a && x < b
            p(i,5) = (2*(x-a))/((c-a)*(b-a));
        elseif x == b
            p(i,5) = 2/(c-a);
        else
            p(i,5) = 2*(c-x)/((c-a)*(c-b));
        end
            
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
        
    end

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

    prob = [ p_A, p_B; p_C, p_D];

%% Importance Sampling Analysis
success = survived_inside;

% f is the indicator function of a successful payload
f = zeros(N,1);
f(success) = 1;

% compute the importance sampling estimate as sum(f.*p./q) / N
E = sum(f.*p./q)/N;


