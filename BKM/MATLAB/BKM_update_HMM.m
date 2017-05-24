function [N, theta, accept] = BKM_update_HMM(N, theta, prior, delta, y, m, f, stdT, N_max)
% A normal random-walk algorithm for the demographic regression parameters,
% Gibbs update for sigma2y

    T = size(N,2);
    %   N1 = N(1,:);
    %   Na = N(2,:);

    D = size(theta,2); 

    fn_BKM_cov = @(xx) BKM_covariates(xx,f,stdT);
    [N_new, accept] = BKM_updateN_HMM(fn_BKM_cov, N, theta, y, delta.N, prior.N, N_max);
    N = N_new;  
    oldlikhood = BKM_calclikhood_HMM(N, theta, y, m, f, stdT, prior.N, N_max);

    % Cycle through each theta (except sigma2) in turn 
    % and propose to update using random walk MH with uniform proposal density:
    
    % For survival/reproduction parameters
    for ii = [1:3,5:7]
        % Keep a record of the current thetaeter value being updated
        oldtheta = theta(ii);
        
        % Propose a new value for the logistic regression parameter theta using a RW with uniform proposal density
        % theta(ii) = runif(1, theta(ii)-delta(ii), theta(ii)+delta(ii));
        theta(ii) = theta(ii) + delta.T(ii)*randn; 
                 
        % Calculate the log(acceptance probability):
        % Calculate the new likelihood value for the proposed move:
        % Calculate the numerator (num) and denominator (den) in turn:
        newlikhood = BKM_calclikhood_HMM(N, theta, y, m, f, stdT, prior.N, N_max);

        % For regression coefficients add in prior (Normal) terms to the acceptance probability
        num = newlikhood - 0.5*(((theta(ii)-prior.T_mu(ii))^2)/prior.T_sigma2(ii));
        den = oldlikhood - 0.5*(((oldtheta-prior.T_mu(ii))^2)/prior.T_sigma2(ii));
        
        % All other prior terms (for other thetas) cancel in the acceptance probability.
        % Proposal terms cancel since proposal distribution is symmetric.

        % Acceptance probability of MH step:
        A = min(1,exp(num-den));

        % To do the accept/reject step of the algorithm        
        % Accept the move with probability A:
        if (rand <= A)  % Accept the proposed move:
            % Update the log(likelihood) value:
            oldlikhood = newlikhood;     
            accept = accept+1;
        else  % Reject proposed move:
            % thetaeter stays at current value:
            theta(ii) = oldtheta;
        end
    end  
    
    % For recovery parameters
    oldlikhood = BKM_recovery(theta, m, f(1:(T-1)), stdT(1:(T-1)));

    for ii = [4,8]
        % Keep a record of the current thetaeter value being updated
        oldtheta = theta(ii);
        
        % Propose a new value for the logistic regression parameter theta using a RW with uniform proposal density
        % theta(ii) = runif(1, theta(ii)-delta(ii), theta(ii)+delta(ii));
        theta(ii) = theta(ii) + delta.T(ii)*randn; 
                 
        % Calculate the log(acceptance probability):
        % Calculate the new likelihood value for the proposed move:
        % Calculate the numerator (num) and denominator (den) in turn:
        newlikhood = BKM_recovery(theta, m, f(1:(T-1)), stdT(1:(T-1)));

        % For regression coefficients add in prior (Normal) terms to the acceptance probability
        num = newlikhood - 0.5*(((theta(ii)-prior.T_mu(ii))^2)/prior.T_sigma2(ii));
        den = oldlikhood - 0.5*(((oldtheta-prior.T_mu(ii))^2)/prior.T_sigma2(ii));
        
        % All other prior terms (for other thetas) cancel in the acceptance probability.
        % Proposal terms cancel since proposal distribution is symmetric.

        % Acceptance probability of MH step:
        A = min(1,exp(num-den));

        % To do the accept/reject step of the algorithm        
        % Accept the move with probability A:
        if (rand <= A)  % Accept the proposed move:
            % Update the log(likelihood) value:
            oldlikhood = newlikhood;     
            accept = accept+1;
        else  % Reject proposed move:
            % thetaeter stays at current value:
            theta(ii) = oldtheta;
        end
    end  
    
    %% For Gibbs update for sigma2
    theta(D) = 1/gamrnd(prior.S(1)+(T-2)/2,1/(prior.S(2)+0.5*sum((y(3:T)-N(3:T)).^2)));  % priorS = [0.001,0.001]
end