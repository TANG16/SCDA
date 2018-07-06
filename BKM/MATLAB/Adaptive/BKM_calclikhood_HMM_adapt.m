function [loglik_comb, loglik_ss, loglik_m] = BKM_calclikhood_HMM_adapt(Na, theta, y, m, f, stdT, priorN, mid, logfact)
    T = length(y);
    loglik_ss = BKM_statespace_HMM_adapt(Na, theta, y, f, stdT, priorN, mid, logfact);
    loglik_m = BKM_recovery(theta, m, f(1:(T-1)), stdT(1:(T-1)));
    loglik_comb = loglik_ss + loglik_m;
end