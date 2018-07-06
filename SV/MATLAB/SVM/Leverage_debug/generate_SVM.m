function [y,h, y0] = generate_SVM(theta,T)
    
    mu = theta(1);
    phi = theta(2);
    sigma2 = theta(3);
    beta = theta(4);
%     rho = (nu-2)/nu;

    h = zeros(1,T);
    h(1,1) = mu + sqrt(sigma2/(1-phi^2))*randn;
    
    sigma = sqrt(sigma2);
    
    for t = 2:T
        h(1,t) = mu + phi*(h(1,t-1) - mu) + sigma*randn;
    end
    
    y0 = exp(h/2).*randn(1,T);
    y = y0 + beta*exp(h);
end