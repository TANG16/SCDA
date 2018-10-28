function loglik = BKM_loglik_N_HMM_adapt_EndCond(Na, Na_prev, Na_prev2, phia, phi1, rho, mid, logfact)

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
              
    loglik = log(sum(P)); % piecewise multiplication enough here
end