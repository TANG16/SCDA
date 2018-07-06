function loglik = loglik_y(y, h, theta)

%     T = length(h);
%     M = size(theta,1);
    
%     mu = theta(:,1);
%     phi = theta(:,2);
%     sigma2 = theta(:,3);
    beta = theta(:,4);

    loglik = -0.5*sum(log(2*pi) + h + ((y-beta*exp(h)).^2)./exp(h));
end