function loglik_ss = BKM_statespace_HMM_NEW(Na, theta, y, f, stdT, priorN, ...
    N_max, logfact)
    T = length(y);
    [phi1, phia, rho, ~] = BKM_covariates(theta,f,stdT);
    sigy = theta(9);
    
    loglam = log(Na) + log(rho) + log(phi1);
    
%     logfact = @(xx) sum(log(1:1:xx));
% %     logfact = @(xx) log(factorial(xx));  % logfact is the log of the factorial: log(x!)
    IND = (0:N_max)';
    loglik = sum(log(binopdf(Na(1:2),priorN(2),priorN(3))));
%     sum(log(poisspdf(N1(1:2),priorN(1)))) 
% 	G = zeros(N_max+1, 1);
%     G(1:N_max,1) = exp(-exp(log(priorN(1))) + IND(1:N_max)*log(priorN(1)) - ...
%                 logfact(IND(1:N_max) + 1)); 
%     G(N_max+1,1) = max(0,1 - sum(G(1:(N_max))));  
%     loglik = loglik + log(sum(G));         

    
    for t = 3:T
        G = zeros(N_max+1, 1);
        P = zeros(N_max+1, 1);
        
        if (t == 3)
            G(1:N_max,1) = exp(-exp(log(priorN(1))) + IND(1:N_max)*log(priorN(1)) - ...
                        logfact(IND(1:N_max) + 1)); 
            G(N_max+1,1) = max(0,1 - sum(G(1:(N_max))));            
        else            
            G(1:N_max,1) = exp(-exp(loglam(t-2)) + IND(1:N_max)*loglam(t-2) - ...
                        logfact(IND(1:N_max) + 1)); 
            G(N_max+1,1) = max(0,1 - sum(G(1:(N_max))));
        end
        %t=4 ==> loglam(2)=log(Na(2)*phi1(2)*rho(2))
        IND_ok = IND((IND + Na(t-1) - Na(t)) > 0);%P(N1(t-1)=k|Na(t-2))
%         if ~isempty(IND_ok)
            P(IND_ok+1) = exp(Na(t)*log(phia(t-1)) + (IND_ok + Na(t-1) - Na(t))*log(1-phia(t-1)) + ...
                    logfact(IND_ok + Na(t-1) + 1) - ...
                    logfact(IND_ok + Na(t-1) - Na(t) + 1) - ...
                    logfact(Na(t) + 1));
%         end
            
        loglik = loglik + log(sum(G.*P)); % piecewise multiplication enough here
    end
  
%     P = ones(N_max+1,1)/(N_max+1);
	G = zeros(N_max+1, 1);
    G(1:N_max,1) = exp(-exp(loglam(T-1)) + IND(1:N_max)*loglam(T-1) - ...
                        logfact(IND(1:N_max) + 1)); 
    G(N_max+1,1) = max(0,1 - sum(G(1:(N_max))));
%     loglik = loglik + log(sum(G .* P));         
    loglik = loglik + log(sum(G));         
    
%     G(1:N_max,1) = exp(-exp(loglam(T)) + IND(1:N_max)*loglam(T) - ...
%                         logfact(IND(1:N_max) + 1)); 
%     G(N_max+1,1) = max(0,1 - sum(G(1:(N_max))));
%     loglik = loglik + log(sum(G .* P));     
            
    loglik_y = sum(- 0.5*(log(2*pi) + log(sigy) + ((y(3:T)-Na(3:T)).^2)/sigy));  
    loglik_ss = loglik + loglik_y;
end