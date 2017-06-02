function [y,h] = generate_SV(theta,T)
    
    mu = theta(1);
    phi = theta(2);
    sigma2 = theta(3);

    h = zeros(1,T);
    h(1,1) = mu + sqrt(sigma2/(1-phi^2))*randn;
    
    sigma = sqrt(sigma2);
    
    for t = 2:T
        h(1,t) = mu + phi*(h(1,t-1) - mu) + sigma*randn;
    end
    
    y = exp(h/2).*randn(1,T);
end