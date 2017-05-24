function [loglik_comb, loglik_ss, loglik_m] = BKM_calclikhood(N, theta, y, m, f, stdT, priorN)
    T = length(y);
    loglik_ss = BKM_statespace(N, theta, y, f, stdT, priorN);
    loglik_m = BKM_recovery(theta, m, f(1:(T-1)), stdT(1:(T-1)));
    loglik_comb = loglik_ss + loglik_m;
end