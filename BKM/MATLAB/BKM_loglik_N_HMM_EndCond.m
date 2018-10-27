function loglik = BKM_loglik_N_HMM_EndCond(Na, phi1, rho, N_max, logfact)


    loglam = log(Na) + log(rho) + log(phi1);
    G = zeros(N_max+1,1);
%     P = ones(N_max+1,1)/(N_max+1);
    
    IND = (0:N_max)';
    G(1:N_max,1) = exp(-exp(loglam) + IND(1:N_max)*loglam - logfact(IND(1:N_max) + 1)); 
    G(N_max+1,1) = max(0,1 - sum(G(1:(N_max))));
%     loglik = log(sum(G  .* P)); % piecewise multiplication enough here
    loglik = log(sum(G)); % piecewise multiplication enough here
end



