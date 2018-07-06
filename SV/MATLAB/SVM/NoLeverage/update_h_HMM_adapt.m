% function [h, accept, A_sum] = update_h_HMM_adapt(y, h, theta, delta_h, mid)
function [h, accept, A_sum] = update_h_HMM_adapt(y, h, theta, delta_h, mid_inv)
% integrate out the odd h(t)'s and impute the even ones
    T = length(y);
    odd = mod(T,2);    
    T2 = floor(T/2);
    T = 2*T2; % to make T even

%     bin_size = bin_midpoint(2)-bin_midpoint(1);
    
    mu = theta(1);
    phi = theta(2);
    sigma2 = theta(3);
    sigma = sqrt(sigma2);
    beta = theta(4);
    
    accept = 0;
    A_sum = 0;
    
    h0 = mu; %mu/(1-phi); % unconditional mean

    Gauss_const = - 0.5*(log(2*pi) + log(sigma2));
   
    for t = 2:2:T % 1:2:T    
%       bin_mid[i,t] <- qnorm(mid[i],phi*(h[t-1] - mu), sigma2_star) 
%       
%       Q_odd[i,t] <- dnorm(y[2*t-1],0, 1.0/exp(bin_mid[i,t]  + mu)) # observation density for odd observations given odd states
%       G_even[i,t] <- dnorm(h[t], mu + phi*(bin_mid[i,t]), sigma2_star) # transition probability to h_even[t] from h_odd[t]       
                      
        %% RW IMPUTATION
        % Keep a record of the current N1 value being updated
        h_old = h(t);
        % normal  RW for h
        h(t) = h_old + delta_h*randn;
        % Calculate the log(acceptance probability):
        % Calculate the new likelihood value for the proposed move:
        % Calculate the numerator (num) and denominator (den) in turn:
        
        % determine the quantiles
        if (t == 2)
%             bin_midpoint = my_norminv(mid,phi*(h0-mu),sigma);
            bin_midpoint = phi*(h0-mu) + sigma.*mid_inv;
        else
%             bin_midpoint = my_norminv(mid,phi*(h(t-2)-mu),sigma);
            bin_midpoint = phi*(h(t-2)-mu) + sigma.*mid_inv;          
        end
        
mu_bin = (mu + bin_midpoint);
exp_mu_bin = exp(mu + bin_midpoint); 

% CURRENTLY UPDATED H BECOMES A QUANTILE AS WELL
bin_midpoint_next = phi*(h(t)-mu) + sigma.*mid_inv;          
bin_midpoint_next_old = phi*(h_old-mu) + sigma.*mid_inv;

        %% Numerator
        num = -0.5*(log(2*pi) + h(t) + ((y(t) - beta*exp(h(t))).^2)/exp(h(t)));
        % integrate out the previous h 
        loglik_int = Gauss_const - 0.5*(((h(t) - mu - phi*bin_midpoint).^2)/sigma2);
        loglik_int = loglik_int - 0.5*(log(2*pi) + mu_bin + ((y(t-1) - beta*exp_mu_bin).^2)./exp_mu_bin);    
   
        num = num + log(sum(exp(loglik_int)));
 
        % integrate out the next h 
        if (t<T)   % t+2 <= T
            loglik_int = Gauss_const - 0.5*(((h(t+2) - mu - phi*bin_midpoint_next).^2)/sigma2);
            loglik_int = loglik_int - 0.5*(log(2*pi) + (mu + bin_midpoint_next) + ((y(t+1) - beta*exp(mu + bin_midpoint_next)).^2)./exp(mu + bin_midpoint_next));    
            num = num + log(sum(exp(loglik_int)));
        end
        
        %% Denominator
        den = -0.5*(log(2*pi) + h_old + ((y(t) - beta*exp(h_old)).^2)/exp(h_old));
        % integrate out the previous h 
        loglik_int = Gauss_const - 0.5*(((h_old - mu - phi*bin_midpoint).^2)/sigma2);
        loglik_int = loglik_int - 0.5*(log(2*pi) + mu_bin + ((y(t-1) - beta*exp_mu_bin).^2)./exp_mu_bin);    
 
        den = den + log(sum(exp(loglik_int)));

        % integrate out the next h 
        if (t<T)   % t+2 <= T
            loglik_int = Gauss_const - 0.5*(((h(t+2) - mu - phi*bin_midpoint_next_old).^2)/sigma2);
            loglik_int = loglik_int - 0.5*(log(2*pi) + (mu + bin_midpoint_next_old) + ((y(t+1) - beta*exp(mu + bin_midpoint_next_old)).^2)./exp(mu + bin_midpoint_next_old));    
            den = den + log(sum(exp(loglik_int)));
        end          
        %% Acceptance Rate
        % Proposal terms cancel since proposal distribution is symmetric.
        % All other prior terms cancel in the acceptance probability. 
        % Acceptance probability of MH step:
        A = min(1,exp(num-den));
 
        A_sum = A_sum + A;
        % To do the accept/reject step of the algorithm:        
        % Accept the move with probability A:
        if (rand <= A)  % Accept the proposed move:
            % Update the log(likelihood) value:
%             oldlikhood = newlikhood;     
            accept = accept+1;
        else  % Reject proposed move:
            % Na stays at current value:
            h(t) = h_old;
        end              
    end
end
