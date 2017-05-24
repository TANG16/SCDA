function [N, theta, newlikhood, accept] = BKM_updateparam(N, theta, oldlikhood, prior, delta, y, m, f, stdT)
% Algorithm 2. Single-update Metropolis–Hastings algorithm using the system process
% (i.e. binomial or Poisson distribution) as the proposal distribution
% A uniform random-walk algorithm for the demographic regression parameters, ?, ?, and ?, and the
% population sizes, N1,0 and Na,0. Alternatively, for ?2 y, we use a Gibbs
% update (conjugate)

    T = size(N,2);
    %   N1 = N(1,:);
    %   Na = N(2,:);

    D = size(theta,2);
    % Uniform RW for the initial states (t=1,2)
    deltaN = delta.N;
    % Uniform RW for regression parameters
    deltaT = delta.T;
    % Conjugate prior on sigma2 (inverse gamma)
    priorS = prior.S; % priorS = [0.001,0.001]
    % Poisson/Binomial priors on initial states (t=1,2) (lam for Poisson for N1; N,p for Bionmial for Na)
    priorN = prior.N; % priorN = [200 2000 0.5]
    % Normal priors (mean and variance) on regression coefficients
    priorT_mu = prior.T_mu; % priorT_mu = 0*ones(D-1,1);
%     priorT_sigma = sqrt(prior.T_sigma2); % prior.T_sigma2 = 100*ones(D-1,1);
    
    fn_BKM_cov = @(xx) BKM_covariates(xx,f,stdT);
    fn_BKM_calclikhood = @(xx) BKM_calclikhood(xx, theta, y, m, f, stdT, priorN);
%     [N_new, newlikhood, accept] = BKM_updateN(fn_BKM_calclikhood, fn_BKM_cov, N, oldlikhood, theta, deltaN, priorN);
%     [N_new, newlikhood, accept] = BKM_updateN_v2(fn_BKM_calclikhood, fn_BKM_cov, N, oldlikhood, theta, deltaN, priorN);
%     [N_new, newlikhood, accept] = BKM_updateN_v3(fn_BKM_calclikhood, fn_BKM_cov, N, oldlikhood, theta, deltaN, priorN);
    [N_new, newlikhood, accept] = BKM_updateN_v4(fn_BKM_calclikhood, fn_BKM_cov, N, oldlikhood, theta, deltaN, priorN);
    N = N_new;
    oldlikhood = newlikhood;
    % Cycle through each theta (except sigma2) in turn 
    % and propose to update using random walk MH with uniform proposal density:
    
    for ii = 1:(D-1)
        % Keep a record of the current thetaeter value being updated
        oldtheta = theta(ii);

        % Propose a new value for the logistic regression parameter theta using a RW with uniform proposal density
        % theta(ii) = runif(1, theta(ii)-delta(ii), theta(ii)+delta(ii));
        theta(ii) = theta(ii) + deltaT(ii)*(2*rand-1); 
                 
        % Calculate the log(acceptance probability):
        % Calculate the new likelihood value for the proposed move:
        % Calculate the numerator (num) and denominator (den) in turn:
        newlikhood = BKM_calclikhood(N, theta, y, m, f, stdT, priorN);

        % For regression coefficients add in prior (Normal) terms to the acceptance probability
%         num = newlikhood + log(normpdf(theta(ii),priorT_mu(ii),priorT_sigma(ii)));
%         den = oldlikhood + log(normpdf(oldtheta,priorT_mu(ii),priorT_sigma(ii)));
%         num = newlikhood - 0.5*(log(2*pi) + log(prior.T_sigma2(ii)) + ((theta(ii)-priorT_mu(ii))^2)/prior.T_sigma2(ii));
%         den = oldlikhood - 0.5*(log(2*pi) + log(prior.T_sigma2(ii)) + ((oldtheta-priorT_mu(ii))^2)/prior.T_sigma2(ii));
        num = newlikhood - 0.5*(((theta(ii)-priorT_mu(ii))^2)/prior.T_sigma2(ii));
        den = oldlikhood - 0.5*(((oldtheta-priorT_mu(ii))^2)/prior.T_sigma2(ii));
        
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
    
    %% Gibbs update for sigma2
    theta(D) = 1/gamrnd(priorS(1)+T/2,1/(priorS(2)+0.5*sum((y-N(2,:)).^2)));  % priorS = [0.001,0.001]
    newlikhood = BKM_calclikhood(N, theta, y, m, f, stdT, priorN);
end