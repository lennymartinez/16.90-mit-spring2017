 function x = logrand(mu, sigma)
% draw sample from a log normal distribution ln N(mu, sigma^2)
z = rand;
x = exp(mu + sigma*z);

end