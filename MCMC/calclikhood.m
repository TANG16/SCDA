function likhood = calclikhood(ni, nj, data, param, nparam, cov, ncov)
%%%%%%
% This function calculates the log likelihood of capture-recapture data
%%%%%%
    % Output: the log(likelihood) value
  
    % First we set the survival and recapture probs:
    % Set up the size of the array containing the survival and recapture probs:
    % Set up the set of cell probabilities (q), and initially set them to be all equal to zero:

%     phi = zeros(1,nj);
%     p = zeros(1,nj);
    q = zeros(ni,nj+1);

    % Set the recapture and survival probs for each year of the study:
    p = param(1)*ones(1,nj);
    exprn = param(2) + param(3)*cov;
    phi = 1./(1+exp(-exprn));
    
%     for (ii = 1:nj) 
%         exprn = param(2) + param(3)*cov(ii);
%         phi(ii) = 1/(1+exp(-exprn));
%         p(ii) = param(1);
%     end

    % Set initial value of the log(likelihood) to be zero
    likhood = 0;

    % Calculate the Multinomial cell probabilities and corresponding contribution to the
    % log(likelihood):
    % Cycle through the number of years that there are releases:
    for ii = 1:ni
        % For diagonal elements:
        q(ii,ii) = phi(ii)*p(ii);
        likhood = likhood + data(ii,ii)*log(q(ii,ii));

        % Calculate the elements above the diagonal:
        if (ii <= (nj-1)) 
            for jj = (ii+1):nj 
                q(ii,jj) = prod(phi(ii:jj))*prod(1-p(ii:(jj-1)))*p(jj);
                likhood = likhood + data(ii,jj)*log(q(ii,jj));
            end
        end
        % Probability of an animal never being seen again
        q(ii,nj+1) = 1 - sum(q(ii,ii:nj));
        likhood = likhood + data(ii,nj+1)*log(q(ii,nj+1));
    end
end