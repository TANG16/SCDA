function [loglik_C, loglik_m, loglik_KF, N_smooth, PHI1, PHIA, RHO, LAMBDA] = BKM_Bugs_combined(theta, y, f, m)
    %% Parameters
    T = max(size(y));
    [T1 , T2] = size(m);
    T2 = T2-1;
    
    TT = size(f,2);
    time = 1:TT;               
    stdT = (time - mean(time))/std(time);
   
    alpha1 = theta(:,1);
    alphaa = theta(:,2);
    alphar = theta(:,3);
    alphal = theta(:,4);
    
    beta1 = theta(:,5);
    betaa = theta(:,6);
    betar = theta(:,7);
    betal = theta(:,8);
    sigy = theta(:,9);
       
	%% Covariates: define the logistic regression equations
    index = alpha1 + beta1*f;
    PHI1 = exp(index)./(1+exp(index));  
    index = alphaa + betaa*f; 
    PHIA = exp(index)./(1+exp(index));  	      
    index = alphar + betar*stdT;
    RHO = exp(index);
    index = alphal + betal*stdT;
    LAMBDA = exp(index)./(1+exp(index));  
	 
	%% Recovery: calculate the cell probabilities for the recovery table
    ind = 1:TT-1;
    phi1 = PHI1(ind);
    phia = PHIA(ind);
    lambda = LAMBDA(ind);
    
    % 	p = zeros(T1,T1);
    % Calculate the diagonal (dies and gets recovered immediately)
    p = diag(lambda.*(1-phi1));	      
    % Calculate value one above the diagonal
    indx = sub2ind(size(p),1:T1-1,2:T1);
    p(indx) = lambda(2:T1).* phi1(1:T1-1).*(1-phia(2:T1));
    
    % Calculate remaining terms above diagonal
    for t1 = 1:(T1-2)
        for t2 = (t1+2):T2
            lphi = sum(log(phia(t1:t2-2)));
            p(t1,t2) = lambda(t2)*phi1(t1)*(1-phia(t2))*exp(lphi);
        end
    end
  	% Probability of an animal never being seen again
	p = [p, 1 - sum(p,2)];  % p(t1, T2+1) = 1 - sum(p(t1,1:T2));
    
%     logp = log(p);
%     logp(~isfinite(logp)) = 0;
%     loglik_m = sum(sum(m.*logp));
    loglik_m = sum(sum(log(mnpdf(m,p))));
    
% Y = mnpdf(X,PROB) returns the pdf for the multinomial distribution with 
% probabilities PROB, evaluated at each row of X.
% X and PROB are m-by-k matrices or 1-by-k vectors, where
% k is the number of multinomial bins or categories. 
% Each row of PROB must sum to one, 
% and the sample sizes for each observation (rows of X) are given by
% the row sums sum(X,2).
% Y is an m-by-k matrix, and mnpdf computes each row of Y 
% using the corresponding rows of the inputs, or replicates them if needed.


    %% Index: state space by normal approximation
    ind = 3:(3+T-1);
    phi1 = PHI1(ind);
    phia = PHIA(ind);
    rho = RHO(ind);   
    
    KFS_params = BKM_SSM(phi1,phia,rho,sigy);
    [N_smooth, ~, nu, F_inv]  = BKM_KFS_multi(y, KFS_params);
    loglik_KF = - 0.5*(T*log(2*pi) - sum(log(abs(F_inv))) + sum(nu.*F_inv.*nu)); 

    loglik_C = loglik_m + loglik_KF;
end