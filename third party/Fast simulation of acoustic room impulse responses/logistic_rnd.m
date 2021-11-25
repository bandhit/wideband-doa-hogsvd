function xx = logistic_rnd(mm,nn,mu,ss)
%logistic_rnd  Generates samples from the logistic distribution
% 
% X = logistic_rnd(M,N,MU,S)
%
% Returns an M-by-N matrix X containing random numbers distributed
% according to the logistic distribution with location (mean) MU and 
% scale S. The code implemented using the inversion method.

% Release date: November 2009
% Author: Eric A. Lehmann, Perth, Australia (www.eric-lehmann.com)

if length(mm)~=1 || length(nn)~=1 || length(mu)~=1 || length(ss)~=1,
    error('all input parameters must be scalars.');
end

yy = rand(mm,nn);

% inverse CDF of the logisitic distribution with location (mean) mu and scale 
% parameter ss, evaluated at yy
xx = mu + ss * log(yy./(1-yy));
