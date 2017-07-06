function logL = SV_loglik_HMM(y, theta, bins, bin_midpoint)

    T = length(y);
    N_bin = length(bin_midpoint);
    % bins are the demeaned volatilities
    mu = theta(1);
    phi = theta(2);
    sigma2 = theta(3);
    sigma = sqrt(sigma2);
    
    G = zeros(N_bin,N_bin);
    % G(ii,jj) = Pr((h_t-mu)in bins(jj:(jj+1)) | (h_t-1 - mu) = bin_midpoint(jj))
    
    for ii = 1:N_bin
        for jj = 1:N_bin
            z_up = (bins(jj+1)- phi*bin_midpoint(ii))/sigma;
            z_dn = (bins(jj)  - phi*bin_midpoint(ii))/sigma;
            G(ii,jj) = normcdf(z_up) - normcdf(z_dn); 
        end
    end
    G = bsxfun(@times,G,1./sum(G,2)); % scale the rows of Gamma to sum up to 1
    delta = ones(1,N_bin)/inv(eye(N_bin) - G + 1); % the stationary distribution, 
                                                   % satisfies delta(Id - Gamma + U) = 1{Nbin,1} 
%     z = exp(-0.5*(bin_midpoint+mu)')*y;
    stdev_y = exp((mu+bin_midpoint)/2);  

    % Algorithm from Section 3.2 in ZMD book: Scaling the likelihood
    % computation
%     log(L(T)) = sum(log(phi(t-1)G(t)P(t))),t=1,...,T)
%     phi(t) = alpha(t)/sum(alpha(t)) 
%     alpha(t) = delta*P(1)*G(1)*...*G(t)*P(t) = alpha(t-1)*G(t)*P(t)
%     L(T) = alpha(T)*ones(m,1);
    alpha = delta;
    loglik = zeros(1,T);
    for t = 1:T
        P = normpdf(y(t),0,stdev_y); 
        alpha = alpha*bsxfun(@times,G,P');
        sumalpha = sum(alpha);
        loglik(t) = log(sumalpha);
        alpha = alpha/sumalpha;        
    end
    logL = sum(loglik);    
end