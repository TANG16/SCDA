function [loglik_comb, loglik_ss, loglik_m] = BKM_calclikhood_HMM(N, theta, y, m, f, stdT, priorN, N_max)
    T = length(y);
    loglik_ss = BKM_statespace_HMM(N, theta, y, f, stdT, priorN, N_max);
    loglik_m = BKM_recovery(theta, m, f(1:(T-1)), stdT(1:(T-1)));
    loglik_comb = loglik_ss + loglik_m;
end