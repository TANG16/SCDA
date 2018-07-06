function [h_mean, h_s, h_var] = h_seq_stats(h, h_mean_old, h_s, t)
    h_mean = h_mean_old + (h - h_mean_old)/t;     
    h_s = h_s + (h - h_mean).*(h - h_mean_old);
    h_var = h_s/(t-1);
    
    % initialize M1 = x1 and S1 = 0.
    % For subsequent x‘s, use the recurrence formulas
    % Mk = Mk-1+ (xk – Mk-1)/k
    % Sk = Sk-1 + (xk – Mk-1)*(xk – Mk).
    % For 2 ? k ? n, the kth estimate of the variance is s2 = Sk/(k – 1).
end