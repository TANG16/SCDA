function N_ess = ESS(x, f)
    % f - truncation parameter
    [N,D] = size(x);
    
    N_ess = zeros(1,D);
    for d = 1:D
        v = sum(autocorr(x(:,d),f));
        N_ess(1,d) = N/(1+2*v);
    end
    
end