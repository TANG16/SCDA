function [y,h, y0] = generate_SVML(theta,T)
    %theta = [-1.00, 0.98, 0.10^2, -0.20, -0.35]; T = 2000

    mu = theta(1);
    phi = theta(2);
    sigma2 = theta(3);
    sigma = sqrt(sigma2);
    beta = theta(4);
    rho = theta(5);

    V = [1, rho*sigma; rho*sigma, sigma^2];
    CV = chol(V);
    E = randn(2,T);
    E = CV'*E;
    epsilon = E(1,:);
    eta = E(2,:);
    
    h = zeros(1,T);
    h(1,1) = mu + sqrt(sigma2/(1-phi^2))*randn;
    
    
    for t = 2:T
        h(1,t) = mu + phi*(h(1,t-1) - mu) + eta(1,t-1);
    end
    
    y0 = exp(h/2).*epsilon;
    y = y0 + beta*exp(h);
end