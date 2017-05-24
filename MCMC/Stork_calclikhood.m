function loglik = Stork_calclikhood(ni, nj, data, param, cov)
%%%%%%
% This function calculates the log likelihood of capture-recapture data
%%%%%%
 
    % Set the recapture and survival probs for each year of the study:
    p = param(1)*ones(1,nj);
    exprn = param(2) + param(3)*cov;
    phi = 1./(1+exp(-exprn));

    % Calculate the Multinomial cell probabilities and corresponding contribution to the
    % log(likelihood):
    % Cycle through the number of years that there are releases:
    
    % For diagonal elements:
    q = diag(phi.*p);	      
    % Calculate the elements above the diagonal:
    for ii = 1:ni 
        for jj = (ii+1):nj 
            q(ii,jj) = prod(phi(ii:jj))*prod(1-p(ii:(jj-1)))*p(jj);
        end
    end
    
    % Probability of an animal never being seen again
	q = [q, 1 - sum(q,2)];   
   
    logq = log(q);
    logq(~isfinite(logq)) = 0;
    loglik = sum(sum(data.*logq));
end