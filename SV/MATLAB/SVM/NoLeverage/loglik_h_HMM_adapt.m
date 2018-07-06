% function loglik = loglik_h_HMM_adapt(y, h, theta, mid)
function loglik = loglik_h_HMM_adapt(y, h, theta, mid_inv)
% if shift integrate out the odd h(t)'s and impute the even ones
% if not shift integrate out the even h(t)'s and impute the odd ones
    T = length(h);
    odd = mod(T,2);    
    T2 = floor(T/2);
    T = 2*T2; % to make T even
    
%     bin_size = bin_midpoint(2)-bin_midpoint(1);
    
    mu = theta(1);
    phi = theta(2);
    sigma2 = theta(3);
    sigma = sqrt(sigma2);
    beta = theta(4);
%     stdev_y = exp((mu+bin_midpoint)/2);  
   
    h0 = mu;  % unconditional mean
    
    Gauss_const = - 0.5*(log(2*pi) + log(sigma2));
        
    loglik = 0;
    for t = 2:2:T 
%       bin_mid[i,t] <- qnorm(mid[i],phi*(h[t-1] - mu), sigma2_star) 
%       
%       Q_odd[i,t] <- dnorm(y[2*t-1],0, 1.0/exp(bin_mid[i,t]  + mu)) # observation density for odd observations given odd states
%       G_even[i,t] <- dnorm(h[t], mu + phi*(bin_mid[i,t]), sigma2_star) # transition probability to h_even[t] from h_odd[t]       
        
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
% G_even        
        % % f(y_t|h_t)
        % loglik = loglik -0.5*(log(2*pi) + h(t) + (y(t)^2)/exp(h(t)));
        % integrate out the previous h 
        % f(h_t|h_t-1)
        loglik_int = Gauss_const - 0.5*(((h(t) - mu - phi*bin_midpoint).^2)/sigma2);
        
% G_odd is now simply a constant for all values as we are using the quantiles of the distribution 
%         % f(h_t-1|h_t-2)
%     %     loglik_int = loglik_int + log(diff(normcdf((bins-phi*(h_prev-mu))/sigma)));   
%         if (t==2)
% %             loglik_int = loglik_int + log(diff(normcdf((bins-phi*(h0-mu))/sigma)));   
%             loglik_int = loglik_int - 0.5*(log(2*pi) + log(sigma2) + ((bin_midpoint - phi*(h0-mu)).^2)/sigma2);
%         elseif (t>2)
% %               loglik_int = loglik_int + log(diff(normcdf((bins-phi*(h(t-2)-mu))/sigma)));   
%             loglik_int = loglik_int - 0.5*(log(2*pi) + log(sigma2) + ((bin_midpoint - phi*(h(t-2)-mu)).^2)/sigma2);
%         end 
% Q_odd     
% % f(y_t-1|h_t-1)
        loglik_int = loglik_int - 0.5*(log(2*pi) + mu_bin + ((y(t-1) - beta*exp_mu_bin).^2)./exp_mu_bin);    
              
        loglik = loglik + log(sum(exp(loglik_int)));   
    end
%     loglik = loglik + T2*log(bin_midpoint(2)-bin_midpoint(1));
end
