function [ll, rec] = loglik_h_leverage(y, h, theta)

    T = length(h);
    
    mu = theta(:,1);
    phi = theta(:,2);
    sigma2 = theta(:,3);
    sigma = sqrt(sigma2);
    beta = theta(4);
    rho = theta(:,5);
        

    rec =  ((h(1) - mu)^2)*(1-phi.^2)*(1-rho^2);
    for t = 2:T
        m = mu + phi*(h(t-1) - mu) + rho*sigma*(y(t-1)-beta*exp(h(t-1)))/exp(h(t-1)/2);
        rec = rec + (h(t) - m)^2;    
    end
    s2 = sigma2*(1-rho^2);
    ll  = rec/s2  - log(1-phi^2)  - log(1-rho^2) + T*log(s2) + T*log(2*pi);
    ll  = -0.5*ll;    
    
end