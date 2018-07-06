function loglik = loglik_h_leverage_HMM_eff(y, h, theta, bins, bin_midpoint)
% integrate out the odd h(t)'s and impute the even ones
    T = length(h);
    odd = mod(T,2);    
    T2 = floor(T/2);
    T = 2*T2; % to make T even
       
    mu = theta(1);
    phi = theta(2);
    sigma2 = theta(3);
    sigma = sqrt(sigma2);
    beta = theta(4);
    rho = theta(5);
    
    h0 = mu;  % unconditional mean
    
    mu_bin = (mu + bin_midpoint);
    exp_mu_bin = exp(mu + bin_midpoint);

    s0 = sigma2/(1-phi^2);
    s2 = sigma2*(1-rho^2);
%     s1 = sqrt(s2);
    Gauss_const0 = - 0.5*(log(2*pi) + log(s0));
    Gauss_const = - 0.5*(log(2*pi) + log(s2));
    
    loglik = zeros(1,T2); % logliks for the integrals consisiting of 3 terms each
    for t = 2:2:T    
        % % f(y_t|h_t)
        % loglik = loglik -0.5*(log(2*pi) + h(t) + (y(t)^2)/exp(h(t)));
        
        % integrate out the previous h 
        % f(h_t|h_t-1)
%         loglik_int = Gauss_const - 0.5*((h(t) - mu - phi*bin_midpoint).^2)/sigma2;
        m = mu + phi*bin_midpoint + ...
            rho*sigma*(y(t-1)-beta*exp_mu_bin)./exp(mu_bin/2);
        loglik_int = Gauss_const - 0.5*((h(t) - m).^2)/s2;
        
%%%% NEW CORRECTION:     
        % % f(y_t-1|h_t-1)
        % loglik_int = loglik_int - 0.5*(log(2*pi) + mu_bin + (y(t-1)^2)./exp_mu_bin);
        % loglik_int = loglik_int - 0.5*(log(2*pi) + (mu + bin_midpoint) + (y(t-1)^2)./exp(mu + bin_midpoint)); 

        loglik_int = loglik_int - 0.5*(log(2*pi) + (mu_bin) + ...
            ((y(t-1) - beta*exp_mu_bin).^2)./exp_mu_bin);          
%%%%

        % % f(h_t-1|h_t-2)
    %     loglik_int = loglik_int + log(diff(normcdf((bins-phi*(h_prev-mu))/sigma))); 
        if (t==2)
%             loglik_int = loglik_int + log(diff(normcdf((bins-phi*(h0-mu))/sigma)));   
%             loglik_int = loglik_int + Gauss_const - 0.5*((bin_midpoint - phi*(h0-mu)).^2)/sigma2;
%             loglik_int = loglik_int  - 0.5*((bin_midpoint - phi*(h0-mu)).^2)/s2;
            loglik_int = loglik_int + Gauss_const0 - 0.5*((bin_midpoint).^2)/s0;
        elseif (t>2)
%               loglik_int = loglik_int + log(diff(normcdf((bins-phi*(h(t-2)-mu))/sigma)));   
%             loglik_int = loglik_int + Gauss_const - 0.5*((bin_midpoint - phi*(h(t-2)-mu)).^2)/sigma2;
            m = phi*(h(t-2)-mu) + rho*sigma*(y(t-2)-beta*exp(h(t-2)))/exp(h(t-2)/2);
            loglik_int = loglik_int + Gauss_const - 0.5*((bin_midpoint - m).^2)/s2;
        end   

%         loglik_int = loglik_int + Gauss_const;     
        loglik(1,t/2) = log(sum(exp(loglik_int)));   
    end
 end
