function loglik_ss = BKM_statespace_HMM(Na, theta, y, f, stdT, priorN, N_max)
% loglik_y
    T = length(y);
    [phi1, phia, rho, ~] = BKM_covariates(theta,f,stdT);
%     lam = N(2:(T-1)).*rho(2:(T-1)).*phi1(2:(T-1));
    sigy = theta(9);
%     G = zeros(N_max+1, T);
%     P = zeros(N_max+1, T);
%     loglam = zeros(T-1,1);   
    loglam = log(Na) + log(rho) + log(phi1);
    
    logfact = @(xx) sum(log(1:1:xx));
%     logfact = @(xx) log(factorial(xx));  % logfact is the log of the factorial: log(x!)
    IND = 0:N_max;
    loglik = sum(log(binopdf(Na(1:2),priorN(1),priorN(2))));
    
    for t = 3:T
        G = zeros(N_max+1, 1);
        P = zeros(N_max+1, 1);
        
        G(1:N_max,1) = exp(-exp(loglam(t-2)) + IND(1:N_max)*loglam(t-2) - logfact(IND(1:N_max))); 
        G(N_max+1,1) = max(0,1 - sum(G(1:(N_max))));

        IND_ok = IND((IND + Na(t-1) - Na(t)) > 0);
        P(IND_ok+1) =  exp(Na(t)*log(phia(t-1)) + (IND_ok + Na(t-1) - Na(t))*log(1-phia(t-1)) + ...
                    arrayfun(logfact,IND_ok + Na(t-1)) - arrayfun(logfact, IND_ok + Na(t-1) - Na(t)) - logfact(Na(t)));
 
%         for ii = 0:(N_max-1)
%             G(ii+1) = exp(-exp(loglam(t-2)) + ii*loglam(t-2) - logfact(ii)); 
% %             G(ii+1) = poisspdf(ii,exp(loglam(t-2))); 
%             if ((ii + Na(t-1) - Na(t)) > 0)
%                 P(ii+1) = exp(Na(t)*log(phia(t-1)) + (ii + Na(t-1) - Na(t))*log(1-phia(t-1)) + ...
%                     logfact(ii + Na(t-1)) - logfact(ii + Na(t-1) - Na(t)) - logfact(Na(t))); 
% %                 P(ii+1) = binopdf(Na(t),ii + Na(t-1),phia(t-1));                 
%             end
%         end
%         G(N_max) = max(0,1 - sum(G(1:(N_max))));
%         if ((N_max + Na(t-1) - Na(t)) > 0)
%             P(N_max+1) = exp(Na(t)*log(phia(t-1)) + (N_max + Na(t-1) - Na(t))*log(1-phia(t-1)) + ...
%                 logfact(N_max + Na(t-1)) - logfact(N_max + Na(t-1) - Na(t)) - logfact(Na(t))); 
% %             P(ii+1) = binopdf(Na(t),ii + Na(t-1),phia(t-1));             
%         end    
        loglik = loglik + log(sum(G.*P)); % piecewise multiplication enough here
    end
    
    loglik_y = sum(- 0.5*(log(2*pi) + log(sigy) + ((y(3:T)-Na(3:T)).^2)/sigy));  
    loglik_ss = loglik + loglik_y;
end