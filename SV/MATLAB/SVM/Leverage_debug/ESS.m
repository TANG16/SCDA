function [N_ess, IF] = ESS(x, f)
    % f - truncation parameter
    [N,D] = size(x);
    N_ess = zeros(1,D);

    if (f > 0)
        for d = 1:D
            v = sum(autocorr(x(:,d),f));
            IF = (1+2*v);
            N_ess(1,d) = N/IF;
        end
    else
        for d = 1:D
            [v, ~, bounds] = autocorr(x(:,d),N-1);
            L =  min(find((v < bounds(1,1)) & (v > bounds(2,1))));
            IF = (1+2*sum(v(1:L)));
            N_ess(1,d) = N/IF;
        end       
    end
end 