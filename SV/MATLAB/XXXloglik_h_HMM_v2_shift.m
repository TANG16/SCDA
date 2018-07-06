function loglik = loglik_h_HMM_v2_shift(y, h, theta, bins, bin_midpoint, shift)
% if shift integrate out the odd h(t)'s and impute the even ones
% if not shift integrate out the even h(t)'s and impute the odd ones
    T = length(h);
    odd = mod(T,2);    
    T2 = floor(T/2);
    T = 2*T2; % to make T even
        
%     bin_size = bin_midpoint(2) - bin_midpoint(1);
        
    mu = theta(1);
    phi = theta(2);
    sigma2 = theta(3);
    stdev_y = exp((mu+bin_midpoint)/2);  
    
    h0 = mu;  % unconditional mean
    if ~shift
        index = 1:2:T; % integrate out h0 (initial) and impute h1
    else
        index = 2:2:T; % integrate out h1 
    end
    
    loglik = 0;
    for t = index  % 2:2:T    
        % loglik = loglik -0.5*(log(2*pi) + h(t) + (y(t)^2)/exp(h(t)));
        % integrate out the previous h 
        loglik_int = - 0.5*(log(2*pi) + log(sigma2) + ((h(t) - mu - phi*bin_midpoint).^2)/sigma2);
        % if (t>1)
            % loglik_int = loglik_int - 0.5*(log(2*pi) + (mu + bin_midpoint) + (y(t-1)^2)./exp(mu + bin_midpoint));    
        % end
    %     loglik_int = loglik_int + log(diff(normcdf((bins-phi*(h_prev-mu))/sigma)));   
        if (t==2)
            loglik_int = loglik_int - 0.5*(log(2*pi) + log(sigma2) + ((bin_midpoint - phi*(h0-mu)).^2)/sigma2);
        elseif (t>2)
            loglik_int = loglik_int - 0.5*(log(2*pi) + log(sigma2) + ((bin_midpoint - phi*(h(t-2)-mu)).^2)/sigma2);
        end      
%%%% NEW CORRECTION:    
        if (t>1)
            loglik_int = loglik_int - 0.5*(log(2*pi) + (mu + bin_midpoint) + (y(t-1)^2)./exp(mu + bin_midpoint));      
        end
        loglik = loglik + log(sum(exp(loglik_int)));   
    end        
%     loglik = loglik + T2*log(bin_midpoint(2)-bin_midpoint(1));
    
end
