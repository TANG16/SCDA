function LL = Heron_HMM_loglik_statespace(theta, N_max1, N_max3, y, f, X2, X4, X2_prior, X4_prior)

    logfact = @(x) sum(log(1:x));
    T = length(y);
    
    G1 = zeros((N_max1+1),T);
    P2 = zeros((N_max1+1),T);
    G3 = zeros((N_max3+1),T);
    P4 = zeros((N_max3+1),T);
    Q = zeros((N_max3+1),T);

    alpha = theta(1:4);
    beta = theta(5:8);
    alpharho = theta(9);
    tauy = theta(10);

    ind = alpha(1) + beta(1)*f;
    phi1 = exp(ind)./(1+exp(ind));
    ind = alpha(2) + beta(2)*f;
    phi2 = exp(ind)./(1+exp(ind));
    ind = alpha(3) + beta(3)*f;
    phi3 = exp(ind)./(1+exp(ind));
    ind = alpha(4) + beta(4)*f;
    phi4 = exp(ind)./(1+exp(ind));
    ind = alpharho*ones(T-1,1); 
    rho = exp(ind);


    loglik = zeros(1,T);
    loglam1 = zeros(1,T-1); 

    for t = 1:(T-1)
        loglam1(t) = log(X4(t)) + log(rho(t)) + log(phi1(t));
    end

    % prior for first transition/augmented observations/observation probabilities
    for ii = 0:N_max1
        G1(ii+1,1) = 1/N_max1; % diffuse =itialisation
        P2(ii+1,1) = X2_prior(1);
    end

    for ii = 0:N_max3
        G3(ii+1,1) = 1/N_max3; % diffuse =itialisation
        P4(ii+1,1) = X4_prior(1);
        Q(ii+1,1) = exp(-tauy*((y(1) - (X2(1) + ii + X4(1))).^2)/2)*sqrt(tauy)/sqrt(2*pi);
    end    

    for t = 2:T
        for ii = 0:(N_max1-1)  % X1 (depends only on [imputed] X4)
            G1(ii+1,t)= exp(-exp(loglam1(t-1)) + ii*loglam1(t-1) - logfact(ii));
        end 
        G1((N_max1+1),t) = max(0,1- sum(G1((1:N_max1),t)));

        for ii = 0:(N_max3) % X3 (depends only on [imputed] X2) & y
            if ((X2(t-1) - ii)>0)
                G3(ii+1,t) = exp(ii*log(phi3(t-1)) + (X2(t-1) - ii)*log(1-phi3(t-1)) + logfact(X2(t-1)) - logfact(abs(X2(t-1)-ii)) - logfact(ii));
            else
                G3(ii+1,t) = 0;
            end
            Q(ii+1,t) = exp(-tauy*((y(t) - (X2(t) + ii + X4(t))).^2)/2)*sqrt(tauy)/sqrt(2*pi);
        end 

        for ii = 0:N_max1 % X2 (depends only on [=tegrated] X1)
            if ((ii - X2(t)) > 0)
                P2(ii+1,t) = exp(X2(t)*log(phi2(t-1)) + (ii - X2(t))*log(1-phi2(t-1)) + logfact(ii) - logfact(abs(ii - X2(t))) - logfact(X2(t)));
            else
                P2(ii+1,t) = 0;
            end
        end

        for ii = 0:(N_max3)  % X4 (depends on [=tegrated] X3)
            if ((ii + X4(t-1) - X4(t)) > 0)
                P4(ii+1,t) = exp(X4(t)*log(phi4(t-1)) + (ii + X4(t-1) - X4(t))*log(1-phi4(t-1)) + logfact(ii + X4(t-1)) - logfact(abs(ii + X4(t-1) - X4(t))) - logfact(X4(t)));
            else
                P4(ii+1,t)= 0;
            end       
        end
    end

    sumsum = zeros(1,T);
    for t = 1:T
        sumsum(t) = exp(log(sum(G1(:,t) .* P2(:,t))) + log(sum(G3(:,t) .* P4(:,t) .* Q(:,t))));
        loglik(t) = log(sumsum(t));
    %     PHI(t) = -loglik(t) + C;
    %     dummy(t) = poisspdf(0,PHI(t));
    end

    LL = sum(loglik) ;
end