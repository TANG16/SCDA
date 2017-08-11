function [theta, A] = update_theta_cond(h, theta, hyper)
% update SV model parameters using the full conditional distributions
% following Kim, Shephard, Chib (1999), RES
    T = length(h);
    [~, rec] = loglik_h(h, theta);
    
%     mu = theta(:,1);     % prior: mu ~ normpdf(c, 0, 10);
%     phi = theta(:,2);    % prior: (phi+1)/2 ~ betapdf((phi+1)/2, 20, 1.5);
%     sigma2 = theta(:,3); % prior: 1/sigma2 ~ gampdf(1./s2, 5/2, 0.05/2);

%     prior_const = [-0.5*(log(2*pi)+log(10)),  - log(beta(20, 1.5)),  2.5*log(0.025), -log(gamma(2.5))];
%     logpdf_norm = @(xx) prior_const(1,1) -0.5*(xx.^2)/10;
%     logpdf_beta = @(xx) prior_const(1,2)  + (20-1)*log(x) + (1.5-1)*log(1-xx); 
%     logpdf_invgamma = @(xx) prior_const(1,3) + prior_const(1,4) - (2.5+1)*log(xx) - 0.025./xx;

    %% Gibbs update for sigma2
    theta(3) = 1/gamrnd(hyper.S(1) + T/2, 1/(hyper.S(2)+0.5*rec));  % hyper.S = [5/2, 0.05/2]
    
    %% MH step for phi
    phi_old = theta(2);
    
    temp = sum(((h(1:T-1)-theta(1))).^2);
    mu_phi = sum((h(2:T)-theta(1)).*(h(1:T-1)-theta(1)))/temp; % phi_hat
    sigma_phi = theta(3)/temp;
    theta(2) = mu_phi + sqrt(sigma_phi)*randn;
    
    if (abs(theta(2)) < 1) 
		log_prior_phi = @(xx) (hyper.P(1)-1)*log((1+xx)/2) + (hyper.P(2)-1)*log((1-xx)/2);
        % hyper.P = [20, 3/2]
        num = log_prior_phi(theta(2)) + 0.5*log(1-theta(2)^2) - 0.5*((h(1)-theta(1))^2*(1-theta(2)^2))/theta(3);
        den = log_prior_phi(phi_old) + 0.5*log(1-phi_old^2) - 0.5*((h(1)-theta(1))^2*(1-phi_old^2))/theta(3);
		
        % Acceptance probability of MH step:
        A = min(1,exp(num-den));
    else
        A = 0;
    end
    % To do the accept/reject step of the algorithm        
    % Accept the move with probability A:
    if (rand <= A)  % Accept the proposed move: 
%         accept = 1;
    else % Reject proposed move:
        % phi stays at current value:
        theta(2) = phi_old;
%         accept = 0;
    end
        
        
    %% Gibbs update for mu
% ADD HYPER ON MU?    
%     sigma_up = theta(3)/((T-1)*(1-theta(2))^2 + (1-theta(2)^2));
    sigma_up = ((T-1)*(1-theta(2))^2 + (1-theta(2)^2))/theta(3) + 1/hyper.M;
    sigma_up = 1/sigma_up;
    mu_up = 0;
    for t = 2:T
        mu_up = mu_up + (h(t) - theta(2)*h(t-1));
    end
    mu_up = mu_up*(1-theta(2)) + (1-theta(2)^2);
    mu_up = mu_up*sigma_up/theta(3);
    theta(1) = mu_up + sqrt(sigma_up)*randn;
end