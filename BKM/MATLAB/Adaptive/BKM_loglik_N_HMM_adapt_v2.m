function loglik = BKM_loglik_N_HMM_adapt_v2(Na, Na_prev, Na_prev2, ...
    phia, phi1, rho, mid, logfact)
% newloglik = BKM_loglik_N_HMM(N(t), N(t-1), phi1(t-1), phia(t-1), rho(t-1), N_bin);
%   y = y(t)
%   Na = N(1,t);
%   Na_prev = N(1,t-1);
%   Na_prev2 = N(1,t-2);
%   phia = phia(t-1) 
%   phi1 = phi1(t-2)
%   rho = rho(t-2)
%     loglam(t-1) = log(Na(t-1-1)) + log(rho(t-1-1)) + log(phi1(t-1-1));
%     loglam(t) = log(Na(t-1)) + log(rho(t-1)) + log(phi1(t-1));
% loglam = loglam(t)

%     logfact = @(xx) sum(log(1:1:xx));
% %     logfact = @(xx) log(factorial(xx));

    loglam = log(Na_prev2) + log(rho) + log(phi1);
    midbin = round(exp(loglam) + exp(loglam/2)*mid);    
    N_q = length(mid);
    
    P = zeros(N_q,1);
%     G  = exp(-exp(loglam) + bin*loglam - logfact(bin + 1)); 
     
    IND_ok = ((midbin +  Na_prev - Na) > 0);
    midbin_ok = midbin(IND_ok);
    P(IND_ok) = exp(Na*log(phia) + (midbin_ok + Na_prev - Na)*log(1-phia) + ...
                  logfact(midbin_ok + Na_prev + 1) - ...
                  logfact(midbin_ok + Na_prev - Na + 1) - ...
                  logfact(Na + 1));
%                 arrayfun(logfact,IND_ok + Na_prev) - ...
%                 arrayfun(logfact, IND_ok + Na_prev - Na) - logfact(Na));
            
%     for ii = 0:(N_bin-1)
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
% %     G(N_bin+1,1) = max(0,1 - sum(G(1:(N_bin))));
%     if ((N_bin + Na_prev - Na) > 0)
%         P(N_bin+1) = exp(Na*log(phia) + (N_bin + Na_prev - Na)*log(1-phi) + ...
%             logfact(N_bin + Na_prev) ...
%             - logfact(N_bin + Na_prev - Na) ...
%             - logfact(Na)); 
% %         P(ii+1) = binopdf(Na,N_bin + Na_prev,phia);                 
%     end    
    loglik = log(sum(P)); % piecewise multiplication enough here
end