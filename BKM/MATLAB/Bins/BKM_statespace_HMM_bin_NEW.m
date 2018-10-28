function loglik_ss = BKM_statespace_HMM_bin_NEW(Na, theta, y, f, stdT, priorN, bin, logfact)
    T = length(y);
    [phi1, phia, rho, ~] = BKM_covariates(theta,f,stdT);
    sigy = theta(9);
   
    loglam = log(Na) + log(rho) + log(phi1);
    
%     logfact = @(xx) sum(log(1:1:xx));
% %     logfact = @(xx) log(factorial(xx));  % logfact is the log of the factorial: log(x!)
    N_bin = length(bin);

    loglik = sum(log(binopdf(Na(1:2),priorN(2),priorN(3))));
    
    for t = 3:T
        P2 = zeros(N_bin, 1);       
        if (t == 3) 
            G2 = exp(-exp(log(priorN(1))) + bin*log(priorN(1)) - logfact(bin + 1));            
        else
            G2 = exp(-exp(loglam(t-2)) + bin*loglam(t-2) - logfact(bin + 1)); 
        end
   
        IND_ok = ((bin + Na(t-1) - Na(t)) > 0);
        bin_ok = bin(IND_ok);
        P2(IND_ok) = exp(Na(t)*log(phia(t-1)) + (bin_ok + Na(t-1) - ...
                    Na(t))*log(1-phia(t-1)) + ...
                    logfact(bin_ok + Na(t-1) + 1) - ...
                    logfact(bin_ok + Na(t-1) - Na(t) + 1) - ...
                    logfact(Na(t) + 1));
 
        loglik = loglik + log(sum(G2.*P2)); % piecewise multiplication enough here
    end
    
    G2 = exp(-exp(loglam(T-1)) + bin*loglam(T-1) - logfact(bin + 1));
    loglik = loglik + log(sum(G2));         
    
    loglik_y = sum(- 0.5*(log(2*pi) + log(sigy) + ((y(3:T)-Na(3:T)).^2)/sigy));  
    loglik_ss = loglik + loglik_y;
end