function x = my_norminv(p,mu,sigma)
%NORMINV Inverse of the normal cumulative distribution function (cdf).
%   X = NORMINV(P,MU,SIGMA) returns the inverse cdf for the normal
%   distribution with mean MU and standard deviation SIGMA, evaluated at
%   the values in P.  The size of X is the common size of the input
%   arguments.  A scalar input functions as a constant matrix of the same
%   size as the other inputs.
%
%   Default values for MU and SIGMA are 0 and 1, respectively.
%
%   [X,XLO,XUP] = NORMINV(P,MU,SIGMA,PCOV,ALPHA) produces confidence bounds
%   for X when the input parameters MU and SIGMA are estimates.  PCOV is a
%   2-by-2 matrix containing the covariance matrix of the estimated parameters.
%   ALPHA has a default value of 0.05, and specifies 100*(1-ALPHA)% confidence
%   bounds.  XLO and XUP are arrays of the same size as X containing the lower
%   and upper confidence bounds.
%
%   See also ERFINV, ERFCINV, NORMCDF, NORMFIT, NORMLIKE, NORMPDF,
%            NORMRND, NORMSTAT.

%   References:
%      [1] Abramowitz, M. and Stegun, I.A. (1964) Handbook of Mathematical
%          Functions, Dover, New York, 1046pp., sections 7.1, 26.2.
%      [2] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical
%          Distributions, 2nd ed., Wiley, 170pp.

%   Copyright 1993-2011 The MathWorks, Inc. 

 

    x0 = -sqrt(2).*erfcinv(2*p);
    x = sigma.*x0 + mu; 
end
