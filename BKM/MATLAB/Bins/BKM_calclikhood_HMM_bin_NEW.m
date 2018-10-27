function [loglik_comb, loglik_ss, loglik_m] = BKM_calclikhood_HMM_bin_NEW(N, theta, y, m, f, stdT, priorN, bin, logfact)
    T = length(y);
    loglik_ss = BKM_statespace_HMM_bin_NEW(N, theta, y, f, stdT, priorN, bin, logfact);
    loglik_m = BKM_recovery(theta, m, f(1:(T-1)), stdT(1:(T-1)));
    loglik_comb = loglik_ss + loglik_m;
end