function [y, tau, h, g] = ucsv_sim(theta, T, hyper)

%     tau0 = 0;
    h0 = theta(1);
    g0 = theta(2);
    omegah = theta(3);
    omegag = theta(4);
    
    phih = 0.9; % to stabilise
    phig = 0.6; % to stabilise
    
%     h = randn(1,T);
%     g = randn(1,T);
    h = zeros(1,T);
    g = zeros(1,T);
    tau = zeros(1,T); 


    h(1) = sqrt(hyper.Vh)*randn;
%     h = cumsum(h);
    g(1) = sqrt(hyper.Vg)*randn;
%     g = cumsum(g);
    
%     std_tau = exp((g0 + omegag*g)/2);
%     tau = std_tau.*randn(1,T);
    tau(1) = sqrt(hyper.Vtau*exp(g0 + omegag*g(1)))*randn;
%     tau = cumsum(tau);
    
    
    for t = 2:T
        h(t) = phih*h(t-1) + randn;
        g(t) = phig*g(t-1) + randn;
        tau(t) = tau(t-1) +  exp((g0 + omegag*g(t))/2)*randn;
    end
        
    y = tau + exp((h0+ omegah*h)/2).*randn(1,T);
end
