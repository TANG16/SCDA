function loglik_ss = BKM_statespace_HMM_bin(Na, theta, y, f, stdT, priorN, bin, logfact)
    T = length(y);
    [phi1, phia, rho, ~] = BKM_covariates(theta,f,stdT);
    sigy = theta(9);
   
    loglam = log(Na) + log(rho) + log(phi1);
    
%     logfact = @(xx) sum(log(1:1:xx));
% %     logfact = @(xx) log(factorial(xx));  % logfact is the log of the factorial: log(x!)
    N_bin = length(bin);

    loglik = sum(log(binopdf(Na(1:2),priorN(1),priorN(2))));
    
    for t = 3:T
        P2 = zeros(N_bin, 1);       
        G2 = exp(-exp(loglam(t-2)) + bin*loglam(t-2) - logfact(bin + 1)); 
   
        IND_ok = ((bin + Na(t-1) - Na(t)) > 0);
        bin_ok = bin(IND_ok);
        P2(IND_ok) = exp(Na(t)*log(phia(t-1)) + (bin_ok + Na(t-1) - ...
                    Na(t))*log(1-phia(t-1)) + ...
                    logfact(bin_ok + Na(t-1) + 1) - ...
                    logfact(bin_ok + Na(t-1) - Na(t) + 1) - ...
                    logfact(Na(t) + 1));
 
        loglik = loglik + log(sum(G2.*P2)); % piecewise multiplication enough here
    end
    
    loglik_y = sum(- 0.5*(log(2*pi) + log(sigy) + ((y(3:T)-Na(3:T)).^2)/sigy));  
    loglik_ss = loglik + loglik_y;
end