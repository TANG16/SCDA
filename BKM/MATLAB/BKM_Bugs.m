function BKM_Bugs(y, f, m, T, T1, T2, rel, stdT)
%  MODEL DESCRIPTION
%***********************************************************************
% inits
tauy = 1;
sigy = 1/tauy;
Na = [1000, 1000, 1092.23, 1100.01, 1234.32, 1460.85, 1570.38, 1819.79,...
    1391.27, 1507.60, 1541.44, 1631.21, 1628.60, 1609.33, 1801.68, 1809.08, 1754.74,...
    1779.48, 1699.13, 1681.39, 1610.46, 1918.45, 1717.07, 1415.69, 1229.02, 1082.02,...
    1096.61, 1045.84, 1137.03, 981.1, 647.67, 992.65, 968.62, 926.83, 952.96, 865.64]';
N1 = 400*ones(T,1);
alpha1 = 1;
alphaa = 2;
alphar = -2;
alphal = -4;
beta1 =-2;
betaa = 0.1;
betar = -0.7;
betal = -0.3;
  
    % Define the priors for the logistic regression parameters
% 	alpha1 = 10*randn; % dnorm(0,0.01)
% 	alphaa = 10*randn; % dnorm(0,0.01)
% 	alphar = 10*randn; % dnorm(0,0.01)
% 	alphal = 10*randn; % dnorm(0,0.01)
% 	beta1 = 10*randn; % dnorm(0,0.01)
% 	betaa = 10*randn; % dnorm(0,0.01)
% 	betar = 10*randn; % dnorm(0,0.01)
% 	betal = 10*randn; % dnorm(0,0.01)
  
	% Define the observation error prior
% 	tauy = gamrnd (0.001,0.001);% dgamma(0.001,0.001)
% 	sigy = 1/tauy;

	% Define the logistic regression equations
    phi1 = zeros(T-1,1);
    phia = zeros(T-1,1);
    rho = zeros(T-1,1);
    lambda = zeros(T-1,1);
    
    for t = 1:(T-1)
        index = alpha1 + beta1*f(t);
		phi1(t) = exp(index)/(1+exp(index)); % logit(phi1(t)) = alpha1 + beta1*f(t) % corresponds to the year 1963
		index = alphaa + betaa*f(t); 
        phia(t) = exp(index)/(1+exp(index)); % logit(phia(t)) = alphaa + betaa*f(t)
		
		% log(rho(t)) = alphar + betar*t % We assume here that t=1
        index = alphar + betar*stdT(t);
		rho(t) = exp(index); % log(rho(t)) = alphar + betar*stdT(t) % We assume here that t=1

		% logit(lambda(t)) = alphal + betal*(t+1)
		index = alphal + betal*stdT(t);
        lambda(t) = exp(index)/(1+exp(index)); % logit(lambda(t)) = alphal + betal*stdT(t)
    end

    r = zeros((T-1)-2,1);
	% Define r(t)
    for t = 3:(T-1)
		r(t-2) = (Na(t+1)+N1(t+1))/(Na(t)+N1(t));
    end

	% Define the initial population priors
%     N1 = zeros(T,1);
%     Na = zeros(T,1);
	
    loglik_N1 = 0;
    loglik_Na = 0;   
    sig1 = 1000^2;
    siga = 1000^2;
    for t = 1:2
        loglik_N1 = loglik_N1  - 0.5*(log(2*pi) + log(sig1) + ((N1(t)-200)^2)/(sig1));
        loglik_Na = loglik_Na  - 0.5*(log(2*pi) + log(siga) + ((Na(t)-1000)^2)/(siga));

% 		N1(t)  = 200 +  1000*randn; %dnorm(200,0.000001)
% 		Na(t)  = 1000 + 1000*randn; %dnorm(1000,0.000001)
    end

	% Define the system process for the census/index data using the Normal approximation

    for t = 3:T
		mean1 = rho(t-1)*phi1(t-1)*Na(t-1);
		meana = phia(t-1)*(N1(t-1)+Na(t-1));
		
% 		tau1(t) = 1/(Na(t-1)*rho(t-1)*phi1(t-1));
% 		taua(t) = 1/((N1(t-1)+Na(t-1))*phia(t-1)*(1-phia(t-1)));
		sig1 = sqrt(Na(t-1)*rho(t-1)*phi1(t-1));
		siga = sqrt((N1(t-1)+Na(t-1))*phia(t-1)*(1-phia(t-1)));
		
% 		N1(t) = mean1 + sig1*randn;
% 		Na(t) = meana + siga*randn;
        loglik_N1 = loglik_N1  - 0.5*(log(2*pi) + log(sig1) + ((N1(t)-mean1)^2)/(sig1));
        loglik_Na = loglik_Na  - 0.5*(log(2*pi) + log(siga) + ((Na(t)-meana)^2)/(siga));
    end
	
	% Define the observation process for the census/index data
	loglik_y = 0;
    for t = 3:T
	    %y(t) ~ dnorm(Na(t),tauy)
        loglik_y = loglik_y - 0.5*(log(2*pi) + log(sigy) + ((y(t)-Na(t))^2)/(sigy));
    end

	% Calculate the no. of birds released each year
	rel = zeros(T1,1);
    for t = 1:T1
	 	rel(t) = sum(m(t,:));
    end
	 


	% Calculate the cell probabilities for the recovery table 
% 	p = zeros(T1,T1);
    % Calculate the diagonal (dies and gets recovered immediately)
    p = diag(lambda.*(1-phi1));	      
%     for t1 = 1 : (T1-1)
% 		p(t1, t1) = lambda(t1)*(1-phi1(t1))
%     end

    % Calculate value one above the diagonal
    indx = sub2ind(size(p),1:T1-1,2:T1);
    p(indx) = lambda(2:T1).* phi1(1:T1-1).*(1-phia(2:T1));

    for t1 = 1 : (T1-1)
		% Calculate the diagonal
% 		p(t1, t1) = lambda(t1)*(1-phi1(t1))
		
		% Calculate value one above the diagonal
        indx = sub2ind(size(p),1:T1-1,2:T1);
        p(indx) = lambda(2:T1).* phi1(1:T1-1).*(1-phia(2:T1));
% 		p(t1, t1+1) = lambda(t1+1)* phi1(t1)*(1-phia(t1+1));

		% Calculate remaining terms above diagonal
        for t2 = (t1+2):T2
            for t = (t1+1):(t2-1)
                lphi(t1, t2, t) = log(phia(t));
            end
            % Probabilities in table
            p(t1,t2) = lambda(t2)*phi1(t1)*(1-phia(t2))*exp(sum(lphi(t1,t2,(t1+1):(t2-1))))
        end
%         for t2 = 1:(t1-1)
% 			% Zero probabilities in lower triangle of table
% 			p(t1, t2) = 0;
%         end
		% Probability of an animal never being seen again
		p(t1, T2+1) = 1 - sum(p(t1,1:T2));	
    end

  	% Probability of an animal never being seen again
	p = [p, 1 - sum(p,2)];  % p(t1, T2+1) = 1 - sum(p(t1,1:T2));
    
	% Final row
% 	p(T1,T1) = lambda(T1)*(1-phi1(T1));	
%     for t = 1:(T1-1)
%         p(T1,t) = 0;
%     end
% 	p(T1,T1+1) = 1 - p(T1,T1);

	% Define the recovery likelihood
    loglik_m = 0;
    for t = 1:T1
        loglik_m = loglik_m + log(mnpdf(m(t,:),p(t,:)));
% 		m(t, 1:(T2+1)) ~ dmulti(p(t,:), rel(t))
    end
end