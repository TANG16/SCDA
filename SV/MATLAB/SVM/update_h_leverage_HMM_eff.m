function [h, accept, A_sum, newloglik] = update_h_leverage_HMM_eff(y, h, theta, delta_h,...
    bins, bin_midpoint, oldloglik)
% integrate out the odd h(t)'s and impute the even ones
    T = length(y);
    odd = mod(T,2);    
    T2 = floor(T/2);
    T = 2*T2; % to make T even
    
    mu = theta(1);
    phi = theta(2);
    sigma2 = theta(3);
    sigma = sqrt(sigma2);
    beta = theta(4);
    rho = theta(5);

    accept = 0;
    A_sum = 0;
      
 	mu_bin = (mu + bin_midpoint);
    exp_mu_bin = exp(mu_bin);

    s0 = sigma2/(1-phi^2);    
    s2 = sigma2*(1-rho^2);
    Gauss_const = - 0.5*(log(2*pi) + log(s2));   
    Gauss_const0 = - 0.5*(log(2*pi) + log(s0));   
    
    newloglik = zeros(1,T2);
    
    for t = 2:2:T % 1:2:T    
        %% RW IMPUTATION
        % Keep a record of the current N1 value being updated
        h_old = h(t);
        % normal  RW for h
        h(t) = h_old + delta_h*randn;
        % Calculate the log(acceptance probability):
        % Calculate the new likelihood value for the proposed move:
        % Calculate the numerator (num) and denominator (den) in turn:

        if (t>2)
            m0 = phi*(h(t-2)-mu) + rho*sigma*(y(t-2)-beta*exp(h(t-2)))./exp(h(t-2)/2);
        end
        LL_int1 = - 0.5*(log(2*pi) + mu_bin + ((y(t-1)-beta*exp_mu_bin).^2)./exp_mu_bin);
        m1 = mu + phi*bin_midpoint + ...
            rho*sigma*(y(t-1)-beta*exp_mu_bin)./exp((mu_bin)/2);
        if (t<T)   % t+2 <= T
            LL_int2 = - 0.5*(log(2*pi) + mu_bin + ((y(t+1)-beta*exp_mu_bin).^2)./exp_mu_bin);
            m2 = mu + phi*bin_midpoint + ...
                rho*sigma*(y(t+1)-beta*exp_mu_bin)./exp((mu_bin)/2);
        end
        
        %% Numerator
        num = -0.5*(log(2*pi) + h(t) + ((y(t)-beta*exp(h(t)))^2)/exp(h(t)));
        % integrate out the previous h 
        loglik_int = Gauss_const - 0.5*((h(t) - m1).^2)/s2;        
        loglik_int = loglik_int + LL_int1;   
        if (t==2)
            loglik_int = loglik_int + Gauss_const0 - 0.5*((bin_midpoint).^2)/s0;
        else
            loglik_int = loglik_int + Gauss_const - 0.5*((bin_midpoint - m0).^2)/s2;
        end    
        NNN = log(sum(exp(loglik_int)));
        num = num + NNN;
        
        % integrate out the next h 
        if (t<T)   % t+2 <= T
            loglik_int = Gauss_const - 0.5*(((h(t+2) - m2).^2)/s2);
            loglik_int = loglik_int + LL_int2;    

            m = phi*(h(t)-mu) + rho*sigma*(y(t)-beta*exp(h(t)))./exp(h(t)/2);
            loglik_int = loglik_int + Gauss_const - 0.5*(((bin_midpoint - m).^2)/s2);
            num = num + log(sum(exp(loglik_int)));
        end    
        
        %% Denominator
        den = -0.5*(log(2*pi) + h_old + ((y(t)-beta*exp(h_old))^2)/exp(h_old));
        % integrate out the previous h 
        loglik_int = Gauss_const - 0.5*(((h_old - m1).^2)/s2);
        loglik_int = loglik_int + LL_int1;    
        if (t==2)
            loglik_int = loglik_int + Gauss_const0 - 0.5*((bin_midpoint).^2)/s0;
        else
            loglik_int = loglik_int + Gauss_const - 0.5*((bin_midpoint - m0).^2)/s2;            
        end    
        DDD = log(sum(exp(loglik_int)));
        den = den + DDD;             
        
        % integrate out the next h 
        if (t<T)   % t+2 <= T
%             loglik_int = Gauss_const - 0.5*(((h(t+2) - mu - phi*bin_midpoint).^2)/sigma2);
%             loglik_int = loglik_int + LL_int2; 
%             loglik_int = loglik_int + Gauss_const - 0.5*(((bin_midpoint - phi*(h_old-mu)).^2)/sigma2);
%             den = den + log(sum(exp(loglik_int)));
             den = den + oldloglik(t/2+1);
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
            newloglik(1,t/2) = NNN;            
        else  % Reject proposed move:
            % Na stays at current value:
            h(t) = h_old;
            % loglik stays at current value:
            newloglik(1,t/2) = DDD;            
        end              
    end
end
