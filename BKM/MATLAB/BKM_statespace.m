function loglik_ss = BKM_statespace(N, theta, y, f, stdT, priorN)
% loglik_y
    N1 = N(1,:);
    Na = N(2,:);
    N_tot = N1 + Na;
    
    T = max(size(y));
   	sigy = theta(9);
    [phi1, phia, rho, ~] = BKM_covariates(theta,f,stdT);
  
    lam = Na(2:(T-1)).*rho(2:(T-1)).*phi1(2:(T-1));
%     fn_logfact = @(xx) sum(log(1:xx));   

    loglik_N1 = sum(log(poisspdf(N1(1:2),priorN(1)))) + sum(log(poisspdf(N1(3:T),lam)));
    loglik_Na = sum(log(binopdf(Na(1:2),priorN(2),priorN(3)))) + sum(log(binopdf(Na(3:T),N_tot(2:(T-1)),phia(2:(T-1)))));
%     loglik_N1 = -2*priorN(1) - sum(lam) ...
%         +  priorN(1)*(N1(1) + N1(2)) + sum(N1(3:T).*lam)...
%         - sum(log(1:N1(1))) - sum(log(1:N1(2))) - sum(arrayfun(fn_logfact,N1));    
%     loglik_Na = sum(Na(3:T).*phia(2:(T-1))) + ...
%         sum((N_tot(2:(T-1))-Na(3:T)).*(1-phia(2:(T-1)))) + ...
%         sum(log(nchoosek(N_tot(2:(T-1)),Na(3:T)))) + ...

    loglik_y = sum(- 0.5*(log(2*pi) + log(sigy) + ((y(3:T)-Na(3:T)).^2)/sigy));
    
    loglik_ss = loglik_N1 + loglik_Na + loglik_y;
end