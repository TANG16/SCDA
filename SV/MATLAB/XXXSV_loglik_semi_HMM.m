function loglik = SV_loglik_semi_HMM(y, h, theta, bins, bin_midpoint, notlast)
    % y = y(1:3), time t-1, t, t+1
    tt = length(y);
    mu = theta(1);
    phi = theta(2);
    sigma2 = theta(3);
    sigma = sqrt(sigma2);
 
%     N_bin = length(bin_midpoint);
% 
%     G1 = zeros(N_bin,1);
%     for ii = 1:N_bin
%         z = (h-mu- phi*bin_midpoint(ii))/sigma;
%         G1(ii,1) = normpdf(z);     
%     end
%     G1 = normpdf((h-mu- phi*bin_midpoint)/sigma);
    G1 = -0.5*(log(2*pi) + log(sigma2) + ((h-mu- phi*bin_midpoint).^2)/sigma2);
    if ~(notlast && (tt==2)) % then it is t=1 and there is onyly y(1) and y(2)
        % observation process
        G1 = G1 - 0.5*(log(2*pi) + log(mu + bin_midpoint) + (y(1)^2)./exp(mu + bin_midpoint));
    end
    G1 = exp(G1);
%     P = normpdf(y,0,exp(h/2));    
    logP = -0.5*(log(2*pi) + log(h) + (y(tt-1)^2)/exp(h));
%     loglik = log(P) + log(sum(G1));
    loglik = logP + log(sum(exp(G1)));
    
    if notlast
% %         G3 = zeros(N_bin,1);
% %         for ii = 1:N_bin
% %             z_up = (bins(ii+1)- phi*(h-mu))/sigma;
% %             z_dn = (bins(ii)  - phi*(h-mu))/sigma;
% %             G3(ii,1) = normcdf(z_up) - normcdf(z_dn);
% %         end
%         z = normcdf((bins - phi*(h-mu))/sigma);
%         G3 = diff(z);
%         G3 = log(G3);
        G3 = -0.5*(log(2*pi) + log(sigma2) + ((bin_midpoint - phi*(h-mu)).^2)/sigma2);
        G3 = G3 - 0.5*(log(2*pi) + log(mu + bin_midpoint) + (y(tt)^2)./exp(mu + bin_midpoint));
        loglik = loglik + log(sum(exp(G3)));
    end
end