    
t = T;

    loglam = log(N(t)) + log(rho(t)) + log(phi1(t));
    G = zeros(N_max+1,1);
    P = ones(N_max+1,1)/N_max+1;
    IND = (0:N_max)';
    G(1:N_max,1) = exp(-exp(loglam) + IND(1:N_max)*loglam - logfact(IND(1:N_max) + 1)); 
    G(N_max+1,1) = max(0,1 - sum(G(1:(N_max))));
    loglik = log(sum(G .* P));
 
    BKM_loglik_N_HMM_EndCond(N(t), phi1(t), rho(t), N_max, logfact)
  
    
    
    loglam = log(N(t-1)) + log(rho(t-1)) + log(phi1(t-1));
    G = zeros(N_max+1,1);
    P = ones(N_max+1,1)/N_max+1;
    IND = (0:N_max)';
    G(1:N_max,1) = exp(-exp(loglam) + IND(1:N_max)*loglam - logfact(IND(1:N_max) + 1)); 
    G(N_max+1,1) = max(0,1 - sum(G(1:(N_max))));
    loglik = log(sum(G .* P));
 
    
    BKM_loglik_N_HMM_EndCond(Na(t-1), phi1(t-1), rho(t-1), N_max, logfact)
    
    
 t = T-1;   
 OK
    BKM_loglik_N_HMM_v2(N(t+1), N(t), N(t-1),...
    phia(t), phi1(t-1), rho(t-1), sigy, N_max, logfact);



    loglam = log(N(t)) + log(rho) + log(phi1);
    G = zeros(N_max+1,1);
    P = ones(N_max+1,1)/N_max+1;
    IND = (0:N_max)';
    G(1:N_max,1) = exp(-exp(loglam) + IND(1:N_max)*loglam - logfact(IND(1:N_max) + 1)); 
    G(N_max+1,1) = max(0,1 - sum(G(1:(N_max))));
    loglik = log(sum(G .* P)); 
    
    
    BKM_loglik_N_HMM_EndCond(Na(t), phi1(t), rho(t), N_max, logfact)