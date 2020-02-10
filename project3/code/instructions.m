
load 'biasing_dist';

% draw N samples from the biasing distribution as a Nx10 matrix
samples = mvnrnd(biasing_dist.mean, biasing_dist.cov, N);

% compute pdf of the biasing distribution q(...) so that q is a Nx1 vector 
% where q(i) is the density evaluated at samples(i,:)
q = mvnpdf(samples, biasing_dist.mean, biasing_dist.cov);


% inside your for loop:

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
    
% f is the indicator function of a successful payload
f = zeros(N,1);
f(success) = 1;

% compute the importance sampling estimate as sum(f.*p./q) / N
