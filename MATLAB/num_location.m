function [num_p] = num_location(sigma, mu, MNM)
% lets assume we have MNSC different places on campus
% probability of each person being in those places is a normal distribution
% function with mu = 0 (for simplicity) and sigma
% this function assigns a PDF to a SC
r(1) = abs(normrnd(mu,sigma.gym));
r(2) = abs(normrnd(mu,sigma.library));
r(3) = abs(normrnd(mu,sigma.dining));
for i = 1:length(r)
    ratio(i) = r(i)/sum(r);
    num_p(i) = MNM*ratio(i);
end
end

