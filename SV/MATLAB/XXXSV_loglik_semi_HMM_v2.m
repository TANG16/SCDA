function loglik = SV_loglik_semi_HMM_v2(Y, H, theta, bins, bin_midpoint)
    % Y = [y(1),y(2),y(3)] time t-1, t, t+1
    tt = length(Y);
    len_H = length(H); % H = [h(0),xh(1),h(2),xh(3),h(4)] with h(1) and h(3) integrated out
    
    mu = theta(1);
    phi = theta(2);
    sigma2 = theta(3);
    sigma = sqrt(sigma2);
 
    % "pure" imputation part, no integration: f(y(2)|h(2))
    loglik =  -0.5*(log(2*pi) + H(2) + (Y(tt-1)^2)/exp(H(2)));
    
    loglik = loglik + integrate(mu, phi, sigma2, Y(1), H(1), H(2), bins, bin_midpoint); % integration wrt h(1)    
    loglik = loglik + integrate(mu, phi, sigma2, Y(tt), H(len_H-1), H(len_H), bins, bin_midpoint); % integration wrt h(3)    
end

function loglik_int = integrate(mu, phi, sigma2, y_curr, h_prev, h_next, bins, bin_midpoint)
%     sigma = sqrt(sigma2);
    loglik_int = 0;

    if ~isnan(h_next)
        loglik_int = loglik_int - 0.5*(log(2*pi) + log(sigma2) + ((h_next - mu - phi*bin_midpoint).^2)/sigma2);
    end
    
    if ~isnan(y_curr)
        loglik_int = loglik_int - 0.5*(log(2*pi) + (mu + bin_midpoint) + (y_curr^2)./exp(mu + bin_midpoint));
    end
    
    if ~isnan(h_prev)
    %     loglik_int = loglik_int + log(diff(normcdf((bins-phi*(h_prev-mu))/sigma)));
        loglik_int = loglik_int -0.5*(log(2*pi) + log(sigma2) + ((bin_midpoint - phi*(h_prev-mu)).^2)/sigma2);
    end
    
    loglik_int = sum(exp(loglik_int));
end