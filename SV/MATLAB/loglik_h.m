function [ll, rec, loglik] = loglik_h(h, theta)

    T = length(h);
%     M = size(theta,1);
    
    mu = theta(:,1);
    phi = theta(:,2);
    sigma2 = theta(:,3);
    
%     a1 = mu;
%     P1 = sigma2./(1-phi.^2);
    
%     ll = -Inf*ones(M,1);
%     ind = ((sigma2 > 0) & (abs(phi) < 1));
        
%     ll = zeros(1,T);
%     ll(1,1) =  ((h(1) - mu)^2)*(1-phi.^2)  - log(1-phi^2);
%     for t = 2:T       
%         ll(1,t) =  (h(t) - mu - phi*(h(t-1) - mu))^2;    
%     end
%     ll  = ll/sigma2 + log(sigma2) + log(2*pi);
%     ll  = -0.5*ll;    
        

    rec =  ((h(1) - mu)^2)*(1-phi.^2) ;
    for t = 2:T
        rec = rec + (h(t) - mu - phi*(h(t-1) - mu))^2;    
    end
    ll  = rec/sigma2  - log(1-phi^2) + T*log(sigma2) + T*log(2*pi);
    ll  = -0.5*ll;
    
%     ll =  ((h(1) - mu)^2)*(1-phi.^2)  - log(1-phi^2);
%     for t = 2:T
%         ll = ll + (h(t) - mu - phi*(h(t-1) - mu))^2;    
%     end
%     ll  = ll/sigma2 + T*log(sigma2) + T*log(2*pi);
%     ll  = -0.5*ll;     
    
end