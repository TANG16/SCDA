function [y,h] = generate_SVt(theta,T)
    
    mu = theta(1);
    phi = theta(2);
    sigma2 = theta(3);
    nu = theta(4);
%     rho = (nu-2)/nu;

    h = zeros(1,T);
    h(1,1) = mu + sqrt(sigma2/(1-phi^2))*randn;
    
    sigma = sqrt(sigma2);
    
    for t = 2:T
        h(1,t) = mu + phi*(h(1,t-1) - mu) + sigma*randn;
    end
    
%     y = exp(rho*h/2).*trnd(nu,1,T);
    y = exp(h/2).*trnd(nu,1,T);
end