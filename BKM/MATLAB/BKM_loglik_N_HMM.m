function loglik = BKM_loglik_N_HMM(y, Na, Na_prev, Na_prev2, phia, phi1, rho, sigy, N_max)
% newloglik = BKM_loglik_N_HMM(N(t), N(t-1), phi1(t-1), phia(t-1), rho(t-1), N_max);
%   y = y(t)
%   Na = N(2,t);
%   Na_prev = N(2,t-1);
%   Na_prev2 = N(2,t-2);
%   phia = phia(t-1) 
%     loglam(t-1) = log(Na(t-1-1)) + log(rho(t-1-1)) + log(phi1(t-1-1));
%     loglam(t) = log(Na(t-1)) + log(rho(t-1)) + log(phi1(t-1));
% loglam = loglam(t)
    logfact = @(xx) sum(log(1:1:xx));
%     logfact = @(xx) log(factorial(xx));

    loglam = log(Na_prev2) + log(rho) + log(phi1);
    G = zeros(N_max+1,1);
    P = zeros(N_max+1,1);
    
    IND = 0:N_max;
    G(1:N_max,1) = exp(-exp(loglam) + IND(1:N_max)*loglam - logfact(IND(1:N_max))); 
    G(N_max+1,1) = max(0,1 - sum(G(1:(N_max))));
    
    IND_ok = IND((IND + Na_prev - Na) > 0);
    P(IND_ok+1) =  exp(Na*log(phia) + (IND_ok + Na_prev - Na)*log(1-phia) + ...
                arrayfun(logfact,IND_ok + Na_prev) - arrayfun(logfact, IND_ok + Na_prev - Na) - logfact(Na));
    
%     for ii = 0:(N_max-1)
%         % logfact is the log of the factorial: log(x!)
% %         G(ii+1) = exp(-exp(loglam) + ii*loglam - logfact(ii)); 
% %         G(ii+1) = poisspdf(ii,exp(loglam)); 
%         if ((ii + Na_prev - Na) > 0)
%             P(ii+1) = exp(Na*log(phia(t-1)) + (ii + Na_prev - Na)*log(1-phia(t-1)) + ...
%                 logfact(ii + Na_prev) ...
%                 - logfact(ii + Na_prev - Na) ...
%                 - logfact(Na));
%         end
%     end
% %     G(N_max+1,1) = max(0,1 - sum(G(1:(N_max))));
%     if ((N_max + Na_prev - Na) > 0)
%         P(N_max+1) = exp(Na*log(phia) + (N_max + Na_prev - Na)*log(1-phi) + ...
%             logfact(N_max + Na_prev) ...
%             - logfact(N_max + Na_prev - Na) ...
%             - logfact(Na)); 
% %         P(ii+1) = binopdf(Na,N_max + Na_prev,phia);                 
%     end    
    loglik = log(sum(G .* P)); % piecewise multiplication enough here
    loglik = loglik - 0.5*(log(2*pi) + log(sigy) + ((y-Na).^2)/sigy);
end