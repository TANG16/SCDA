function [ll, rec] = loglik_h(h, theta)

    T = length(h);
    M = size(theta,1);
    
    mu = theta(:,1);
    phi = theta(:,2);
    sigma2 = theta(:,3);
    
%     a1 = mu;
%     P1 = sigma2./(1-phi.^2);
    
    ll = -Inf*ones(M,1);
    
    ind = ((sigma2 > 0) & (abs(phi) < 1));
        
    rec =  ((h(1) - mu(ind)).^2).*(1-phi(ind).^2);
    for t = 2:T
        mm = (h(t) - mu(ind) - phi(ind).*(h(t-1) - mu(ind)));
        rec = rec + (mm.^2);    
    end
    ll(ind) = rec./sigma2(ind) + T*log(sigma2(ind)) - log(1-phi(ind).^2);
    ll(ind) = ll(ind) + T*log(2*pi);
    ll(ind) = -0.5*ll(ind);    
    
%     ll(ind) = log(P1(ind)) + ((h(1) - a1(ind)).^2)./P1(ind);
%     ll(ind) = ll(ind) + (T-1)*log(sigma2(ind));
%     for t = 2:T
%         mm = (h(t) - mu(ind) - phi(ind).*(h(t-1) - mu(ind)));
%         ll(ind) = ll(ind) + (mm.^2)./sigma2(ind);    
%     end
%     ll(ind) = ll(ind) + T*log(2*pi);
%     ll(ind) = -0.5*ll(ind);
end