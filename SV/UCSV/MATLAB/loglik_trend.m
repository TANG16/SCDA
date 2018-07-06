function ll = loglik_trend(y, tau, h, g, theta, hyper)
 
    h0 = theta(:,1);
    g0 = theta(:,2);
    omegah = theta(:,3);
    omegag = theta(:,4);

    mean_y = tau;
    var_y = exp(h0 + omegah*h);
    mean_tau = [0, tau(1:(end-1))];
    var_tau = exp(g0 + omegag*g);
    var_tau(1) = var_tau(1)*hyper.Vtau;
       
    ll_y = -0.5*(log(2*pi) + log(var_y) + ((y-mean_y).^2)/var_y);
    ll_tau = -0.5*(log(2*pi) + log(var_tau) + ((tau-mean_tau).^2)/var_tau);
    
    ll = sum(ll_y) + sum(ll_tau);
end