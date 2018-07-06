% function loglik = loglik_h_HMM_adapt(y, h, theta, mid)
function loglik = loglik_h_HMM_adapt_eff(y, h, theta, mid_inv)
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
    s0 = sigma2/(1-phi^2);

%     stdev_y = exp((mu+bin_midpoint)/2);  
   
    h0 = mu;  % unconditional mean

    Gauss_const = - 0.5*(log(2*pi) + log(sigma2));
    Gauss_const0 = - 0.5*(log(2*pi) + log(s0));

    loglik = zeros(1,T2); % logliks for the integrals consisiting of 3 terms each
    for t = 2:2:T     
        % determine the quantiles
        if (t == 2)
% %             bin_midpoint = my_norminv(mid,phi*(h0-mu),sigma);
%             bin_midpoint = phi*(h0-mu) + sigma.*mid_inv;
            bin_midpoint = sqrt(s0).*mid_inv;
        else
%             bin_midpoint = my_norminv(mid,phi*(h(t-2)-mu),sigma);
            bin_midpoint = phi*(h(t-2)-mu) + sigma.*mid_inv;          
        end
% G_even        
        % % f(y_t|h_t)
        % loglik = loglik -0.5*(log(2*pi) + h(t) + (y(t)^2)/exp(h(t)));
        % integrate out the previous h 
        % f(h_t|h_t-1)
        loglik_int = - 0.5*(((h(t) - mu - phi*bin_midpoint).^2)/sigma2);
        
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
        loglik_int = loglik_int - 0.5*(log(2*pi) + (mu + bin_midpoint) + (y(t-1)^2)./exp(mu + bin_midpoint));    
        loglik_int = loglik_int + Gauss_const;                 
        loglik(1,t/2) = log(sum(exp(loglik_int)));   
    end
end
