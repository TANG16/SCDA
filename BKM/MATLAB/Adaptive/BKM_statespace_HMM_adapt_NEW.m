function loglik_ss = BKM_statespace_HMM_adapt_NEW(Na, theta, y, f, stdT, priorN, mid, logfact)
% loglik_y
    T = length(y);
    [phi1, phia, rho, ~] = BKM_covariates(theta,f,stdT);
%     lam = N(2:(T-1)).*rho(2:(T-1)).*phi1(2:(T-1));
    sigy = theta(9);
    loglam = log(Na) + log(rho) + log(phi1);
    
%     logfact = @(xx) sum(log(1:1:xx));
% %     logfact = @(xx) log(factorial(xx));  % logfact is the log of the factorial: log(x!)
    N_mid = length(mid);
    
% loglik = zeros(T,1); 
% loglik(1:2) = sum(log(binopdf(Na(1:2),priorN(1),priorN(2))));
    loglik = sum(log(binopdf(Na(1:2),priorN(2),priorN(3))));
    
    for t = 3:T
        
%       midbin[i+1,t] = qnorm(mid[i+1],exp(loglam[t-1]),1/exp(loglam[t-1]))
        if (t == 3)
            midbin = (priorN(1) + sqrt(priorN(1))*mid); % N1(t-1)
        else
            midbin = (exp(loglam(t-2)) + sqrt(exp(loglam(t-2)))*mid); % N1(t-1)
        end
% G is now simply a constant for all values as we are using the quantiles of the distribution 
% and a normal approximation (so they are exact).
         
        P2 = zeros(N_mid, 1);      
%         G2 = exp(-exp(loglam(t-2)) + mid*loglam(t-2) - ...
%                             logfact(bin + 1)); 
 
        IND_ok = ((midbin + Na(t-1) - Na(t)) > 0);
        mid_ok = midbin(IND_ok);
        P2(IND_ok) = exp(Na(t)*log(phia(t-1)) + (round(mid_ok + Na(t-1)) - Na(t))*log(1-phia(t-1)) + ...
                    logfact(round(mid_ok + Na(t-1)) + 1) - ...
                    logfact(round(mid_ok + Na(t-1)) - Na(t) + 1) - ...
                    logfact(Na(t) + 1));
%         loglik(t) = log(sum(P2));
        loglik = loglik + log(sum(P2));
    end
    
%     midbin = (exp(loglam(T-1)) + sqrt(exp(loglam(T-1)))*mid); % N1(T)

    loglik_y = sum(- 0.5*(log(2*pi) + log(sigy) + ((y(3:T)-Na(3:T)).^2)/sigy));  
    loglik_ss = loglik + loglik_y;
end