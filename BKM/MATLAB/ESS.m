function N_ess = ESS(x, f)
    % f - truncation parameter
    [N,D] = size(x);
    N_ess = zeros(1,D);
    
%     if ~exist('autocorr','builtin')
%         autocorr = @(xx,ll) acf(xx,ll);
%     end

    if (f > 0)
        for d = 1:D
            v = sum(autocorr(x(:,d),f));
            N_ess(1,d) = N/(1+2*v);
        end
    else
        for d = 1:D
            [v, ~, bounds] = autocorr(x(:,d),N-1);
            L =  min(find((v < bounds(1,1)) & (v > bounds(2,1))));
            N_ess(1,d) = N/(1+2*sum(v(1:L)));
        end       
    end
end 