function [param, likhood, accept] = Stork_updateparam(nparam, param, cov, ni, nj, data, likhood, alphap, betap, mu, sig, delta)
%%%%%%
% Function for updating the parameters values:
%%%%%%
    
  % output = updateparam(nparam, param, ncov, cov, ni, nj, data, likhood, alphap, betap, mu, sig, delta)
  % nparam = 3
  % ncov = 1
  % delta = [0.05,0.1,0.1)
  % param = [0.9,0.7,0.1)
  % % Beta prior on recapture probability
  % alphap = 1
  % betap = 1
  % % Normal priors (mean and variance) on regression coefficients
  % mu = [0,0)
  % sig2 = [10,10)
  % sig = sqrt(sig2)
    accept = 0; 
    % Cycle through each parameter in turn and propose to update using
    % random walk MH with Uniform proposal density:
    for ii = 1:nparam 
        % Keep a record of the current parameter value being updated
        oldparam = param(ii);

        % Propose a new value for the parameter using a random walk with
        % Uniform proposal density
%         param(ii) = runif(1, param(ii)-delta(ii), param(ii)+delta(ii));
        param(ii) = param(ii) + 2*delta(ii)*rand;

        % Automatically reject any moves where recapture prob is outside [0,1]
        if (param(1) >= 0 && param(1) <= 1) 
            % Calculate the log(acceptance probability):
    
            % Calculate the new likelihood value for the proposed move:
            % Calculate the numerator (num) and denominator (den) in turn:
            newlikhood = Stork_calclikhood(ni, nj, data, param, cov);

            if (ii == 1) 
                % For recapture probability add in prior (Beta) terms to the acceptance probability
                num = newlikhood + log(betapdf(param(ii),alphap,betap));
                den = likhood + log(betapdf(oldparam,alphap,betap));
            else 
                % For regression coefficients add in prior (Normal) terms to the acceptance probability
                num = newlikhood + log(normpdf(param(ii),mu(ii-1),sig(ii-1)));
                den = likhood + log(normpdf(oldparam,mu(ii-1),sig(ii-1)));
            end

            % All other prior terms (for other parameters) cancel in the acceptance probability.
            % Proposal terms cancel since proposal distribution is symmetric.
            
            % Acceptance probability of MH step:
            A = min(1,exp(num-den));
        else 
            A = 0;
        end

        % To do the accept/reject step of the algorithm:
        u = rand; % Simulate a random number in [0,1]:
        
        % Accept the move with probability A:
        if (u <= A)  % Accept the proposed move:
            % Update the log(likelihood) value:
            likhood = newlikhood;     
            accept = accept+1;
        else  % Reject proposed move:
            % parameter stays at current value:
            param(ii) = oldparam;
        end
    end
end