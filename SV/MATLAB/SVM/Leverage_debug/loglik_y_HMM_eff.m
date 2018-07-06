function loglik = loglik_y_HMM_eff(y, h, theta, bins, bin_midpoint)
% if shift integrate out the odd h(t)'s and impute the even ones
% if not shift integrate out the even h(t)'s and impute the odd ones
    T = length(h);
    odd = mod(T,2);    
    T2 = floor(T/2);
    T = 2*T2; % to make T even
    
%     bin_size = bin_midpoint(2)-bin_midpoint(1);
    
%     mu = theta(1);
%     phi = theta(2);
%     sigma2 = theta(3);
    beta = theta(4);
  
%     h0 = mu;  % unconditional mean
%     mu_bin = (mu + bin_midpoint);
%     exp_mu_bin = exp(mu + bin_midpoint);
%     Gauss_const = - 0.5*(log(2*pi) + log(sigma2));
        
%     loglik = 0;

    loglik = sum( - 0.5*(log(2*pi) + h(2:2:T) + ((y(2:2:T)-beta*exp(h(2:2:T))).^2)./exp(h(2:2:T))) );         
    
%     for t = 2:2:T    
%         % integrate out the previous h 
%         % f(h_t|h_t-1)
%         loglik_int = Gauss_const - 0.5*(((h(t) - mu - phi*bin_midpoint).^2)/sigma2);
%         % % f(y_t|h_t)
%         loglik_int = loglik_int - 0.5*(log(2*pi) + h(t) + ((y(t)-beta*exp(h(t)))^2)/exp(h(t))) ;         
%         % f(h_t-1|h_t-2)
%         if (t==2)
%             loglik_int = loglik_int + Gauss_const - 0.5*(((bin_midpoint - phi*(h0-mu)).^2)/sigma2);
%         elseif (t>2)
%             loglik_int = loglik_int + Gauss_const - 0.5*(((bin_midpoint - phi*(h(t-2)-mu)).^2)/sigma2);
%         end 
% %%%% NEW CORRECTION:        
%         % % f(y_t-1|h_t-1)
%         loglik_int = loglik_int - 0.5*(log(2*pi) + mu_bin + ((y(t-1)-beta*exp_mu_bin).^2)./exp_mu_bin) ;         
%         loglik = loglik + log(sum(exp(loglik_int)));   
%     end
end
