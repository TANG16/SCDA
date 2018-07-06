function [N, y] = BKM_genarate
    sc = 1;

    [y_true, T, time, stdT, f, m, T1, T2] = BKM_Data_HMM(sc);

    theta = [0.5237    1.4728   -1.0636   -4.5935 ...
        -0.1528   -0.2435   -0.3597   -0.3546 ...
           6.0304e+04];
    [phi1, phia, rho, ~] = BKM_covariates(theta,f,stdT);

    y = zeros(1,T);
    N = zeros(2,T);
    % Na = zeros(1,T2);

    prior.N = [200/sc 2000/sc 0.5];

    N(1,1:2) = poissrnd(prior.N(1),1,2);
    N(2,1:2) = binornd(prior.N(2),prior.N(3));

    for t = 3:T
        lam = N(2,t-1).*rho(t-1).*phi1(t-1);
        N(1,t) = poissrnd(lam);
        N(2,t) = binornd(N(1,t-1)+N(2,t-1),phia(t-1));
    end 
    y(1,3:T) = N(2,3:T) + sqrt(theta(9))*randn(1,T-2);

end
